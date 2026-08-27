#!/usr/bin/env python3
# merge_gtexv8_mashr_predictdb_eqtl.py
#
# Merge multiple GTEx v8 PredictDB / MetaXcan mashr sqlite DB files.
#
# Key rule:
#   A model is tissue × gene, not gene alone.
#   Therefore both weights.gene and extra.gene are rewritten as:
#       <tissue>.<ENSG...>
#
# Example:
#   python merge_gtexv8_mashr_predictdb.py \
#     --input_dir /path/to/gtex_v8/eqtl \
#     --pattern "mashr_*.db" \
#     --output combine_gtexv8_mashr_eqtl.db

import argparse
import glob
import os
import re
import sqlite3
from pathlib import Path

import pandas as pd


def infer_tissue_from_filename(path: str) -> str:
    """Infer tissue name from GTEx v8 PredictDB file name."""
    name = Path(path).name
    name = re.sub(r"\.db$", "", name)
    name = re.sub(r"^mashr_", "", name)
    name = re.sub(r"^eqtl_", "", name)
    return name


def list_tables(conn):
    q = "SELECT name FROM sqlite_master WHERE type='table'"
    return [x[0] for x in conn.execute(q).fetchall()]


def table_exists(conn, table):
    return table.lower() in [t.lower() for t in list_tables(conn)]


def read_table(conn, table):
    return pd.read_sql_query(f'SELECT * FROM "{table}"', conn)


def write_table(df, conn, table, if_exists="append"):
    df.to_sql(table, conn, if_exists=if_exists, index=False)


def prefix_gene(gene, tissue):
    gene = str(gene)
    # Avoid double-prefixing if the db was already processed.
    if gene.startswith(tissue + "."):
        return gene
    return f"{tissue}.{gene}"


def create_index(conn, table, cols, index_name):
    cols_sql = ", ".join([f'"{c}"' for c in cols])
    conn.execute(f'CREATE INDEX IF NOT EXISTS "{index_name}" ON "{table}" ({cols_sql})')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Folder containing GTEx v8 mashr .db files")
    parser.add_argument("--pattern", default="mashr_*.db", help='Glob pattern, default: "mashr_*.db"')
    parser.add_argument("--output", required=True, help="Output merged sqlite db")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite output db if it exists")
    args = parser.parse_args()

    db_files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not db_files:
        raise SystemExit(f"No db files found: {args.input_dir}/{args.pattern}")

    out = Path(args.output)
    if out.exists():
        if args.overwrite:
            out.unlink()
        else:
            raise SystemExit(f"Output exists: {out}. Use --overwrite or delete it first.")

    out_conn = sqlite3.connect(str(out))

    map_rows = []
    summary_rows = []

    first_weights = True
    first_extra = True
    first_sample_info = True

    for db in db_files:
        tissue = infer_tissue_from_filename(db)
        print(f"[INFO] Processing {tissue}: {db}")

        in_conn = sqlite3.connect(db)
        tables = list_tables(in_conn)

        # ---- weights ----
        if not table_exists(in_conn, "weights"):
            print(f"[WARN] No weights table in {db}; skipping")
            in_conn.close()
            continue

        weights = read_table(in_conn, "weights")
        if "gene" not in weights.columns:
            raise ValueError(f"{db}: weights table has no 'gene' column")

        original_weight_genes = weights["gene"].astype(str).copy()
        weights["gene"] = weights["gene"].map(lambda x: prefix_gene(x, tissue))

        write_table(weights, out_conn, "weights", if_exists="replace" if first_weights else "append")
        first_weights = False

        # Mapping from prefixed model ID back to tissue and Ensembl gene.
        tmp_map = pd.DataFrame({
            "model_id": weights["gene"].astype(str),
            "tissue": tissue,
            "gene": original_weight_genes.astype(str)
        }).drop_duplicates()
        map_rows.append(tmp_map)

        n_weight_rows = len(weights)
        n_weight_models = weights["gene"].nunique()

        # ---- extra ----
        # PredictDB/MetaXcan dbs usually have 'extra' table.
        # SQLite is case-insensitive for table names, but pandas query needs the actual spelling.
        extra_table = None
        for t in tables:
            if t.lower() == "extra":
                extra_table = t
                break

        n_extra_rows = 0
        n_extra_models = 0
        if extra_table is not None:
            extra = read_table(in_conn, extra_table)
            if "gene" in extra.columns:
                extra["gene"] = extra["gene"].map(lambda x: prefix_gene(x, tissue))
            n_extra_rows = len(extra)
            n_extra_models = extra["gene"].nunique() if "gene" in extra.columns else 0

            # Use canonical lowercase table name "extra" in merged DB.
            write_table(extra, out_conn, "extra", if_exists="replace" if first_extra else "append")
            first_extra = False
        else:
            print(f"[WARN] No extra table in {db}")

        # ---- sample_info, optional ----
        if table_exists(in_conn, "sample_info"):
            sample_info = read_table(in_conn, "sample_info")
            sample_info.insert(0, "tissue", tissue)
            write_table(
                sample_info,
                out_conn,
                "sample_info",
                if_exists="replace" if first_sample_info else "append",
            )
            first_sample_info = False

        summary_rows.append({
            "db_file": db,
            "tissue": tissue,
            "n_weight_rows": n_weight_rows,
            "n_weight_models": n_weight_models,
            "n_extra_rows": n_extra_rows,
            "n_extra_models": n_extra_models,
        })

        in_conn.close()

    if map_rows:
        model_map = pd.concat(map_rows, ignore_index=True).drop_duplicates()
        write_table(model_map, out_conn, "gene_tissue_map", if_exists="replace")

    summary = pd.DataFrame(summary_rows)
    write_table(summary, out_conn, "merge_summary", if_exists="replace")

    # Helpful indexes for downstream queries.
    create_index(out_conn, "weights", ["gene"], "idx_weights_gene")
    if table_exists(out_conn, "extra"):
        create_index(out_conn, "extra", ["gene"], "idx_extra_gene")
    create_index(out_conn, "gene_tissue_map", ["model_id"], "idx_gene_tissue_map_model_id")
    create_index(out_conn, "gene_tissue_map", ["gene"], "idx_gene_tissue_map_gene")

    out_conn.commit()

    # Simple consistency checks.
    n_weights_gene = out_conn.execute("SELECT COUNT(DISTINCT gene) FROM weights").fetchone()[0]
    n_extra_gene = out_conn.execute("SELECT COUNT(DISTINCT gene) FROM extra").fetchone()[0] if table_exists(out_conn, "extra") else 0
    n_intersection = out_conn.execute("""
        SELECT COUNT(DISTINCT e.gene)
        FROM extra e
        INNER JOIN weights w ON e.gene = w.gene
    """).fetchone()[0] if table_exists(out_conn, "extra") else 0

    print("")
    print("[DONE] Merged DB:", out)
    print("[CHECK] distinct weights.gene:", n_weights_gene)
    print("[CHECK] distinct extra.gene:", n_extra_gene)
    print("[CHECK] extra ∩ weights:", n_intersection)
    print("[INFO] Added tables: gene_tissue_map, merge_summary")

    out_conn.close()


if __name__ == "__main__":
    main()
