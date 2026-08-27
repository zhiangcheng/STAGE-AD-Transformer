#!/usr/bin/env python3
"""Step 01: generate a synthetic, model-ready SexReg-AD dataset."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.constants import (
    CELL_TYPES,
    CLASS_TO_ID,
    FEMALE_GENES,
    KNOWN_AD_GENES,
    MALE_GENES,
    MECHANISMS,
    TISSUES,
    X_GENES,
)
from sexreg_ad.utils import ensure_dir, neglog10, p_from_z


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic sex-aware AD regulatory features."
    )
    parser.add_argument("--n-variants", type=int, default=300_000)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--seed", type=int, default=20260618)
    return parser.parse_args()


def chromosome_weights() -> tuple[list[str], np.ndarray]:
    chromosomes = [str(index) for index in range(1, 23)] + ["X"]
    weights = np.array(
        [
            249,
            242,
            198,
            190,
            181,
            171,
            159,
            146,
            141,
            135,
            135,
            134,
            115,
            107,
            102,
            90,
            83,
            80,
            59,
            64,
            47,
            51,
            156,
        ],
        dtype=float,
    )
    return chromosomes, weights / weights.sum()


def expected_targets(
    labels: np.ndarray, is_x: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assign synthetic tissue, cell, and mechanism targets by label group."""

    n_rows = len(labels)
    tissue = np.empty(n_rows, dtype=int)
    cell = np.empty(n_rows, dtype=int)
    mechanism = np.empty(n_rows, dtype=int)

    rules = [
        (is_x, [0, 1, 2, 4], [0, 1, 4, 6], [7]),
        (
            (~is_x) & (labels == "female_biased"),
            [0, 1, 2, 4, 5],
            [0, 1, 6],
            [0, 1, 4, 5],
        ),
        (
            (~is_x) & (labels == "male_biased"),
            [6, 7, 8, 0],
            [2, 3, 5],
            [2, 3, 4],
        ),
        (
            (~is_x) & (labels == "sex_interaction"),
            [0, 2, 4, 6],
            [0, 1, 3],
            [0, 2, 5, 7],
        ),
        (
            (~is_x) & (labels == "shared"),
            [0, 1, 2, 4],
            [0, 4],
            [0, 5, 6],
        ),
        (
            (~is_x) & (labels == "null"),
            list(range(len(TISSUES))),
            list(range(len(CELL_TYPES))),
            list(range(len(MECHANISMS))),
        ),
    ]
    for mask, tissue_choices, cell_choices, mechanism_choices in rules:
        count = int(mask.sum())
        tissue[mask] = rng.choice(tissue_choices, size=count)
        cell[mask] = rng.choice(cell_choices, size=count)
        mechanism[mask] = rng.choice(mechanism_choices, size=count)
    return tissue, cell, mechanism


def simulate_chunk(
    start: int, size: int, rng: np.random.Generator
) -> pd.DataFrame:
    chromosomes, chromosome_probabilities = chromosome_weights()
    chromosome = rng.choice(chromosomes, size=size, p=chromosome_probabilities)
    chromosome_lengths = {
        str(index): int(3e9 * chromosome_probabilities[index - 1])
        for index in range(1, 23)
    }
    chromosome_lengths["X"] = 156_000_000
    bp = np.zeros(size, dtype=np.int64)
    for name in np.unique(chromosome):
        mask = chromosome == name
        bp[mask] = rng.integers(1, chromosome_lengths[name], mask.sum())

    variant_id = np.array(
        [f"{chrom}:{position}:A:G" for chrom, position in zip(chromosome, bp)],
        dtype=object,
    )
    rsid = np.array(
        [f"rsSIMV3{start + index:010d}" for index in range(size)], dtype=object
    )
    maf = np.clip(rng.beta(0.8, 4.5, size), 0.001, 0.5)
    info = np.clip(rng.normal(0.965, 0.035, size), 0.70, 1)
    ld_score = np.exp(rng.normal(1.6, 0.75, size))
    distance_to_tss = np.clip(np.exp(rng.normal(4.2, 1.25, size)), 0, 10_000)

    label_draw = rng.random(size)
    labels = np.full(size, "null", dtype=object)
    labels[label_draw < 0.0015] = "shared"
    labels[(label_draw >= 0.0015) & (label_draw < 0.0022)] = "female_biased"
    labels[(label_draw >= 0.0022) & (label_draw < 0.0028)] = "male_biased"
    labels[(label_draw >= 0.0028) & (label_draw < 0.0032)] = "sex_interaction"

    is_x = chromosome == "X"
    x_draw = rng.random(size)
    labels[is_x & (x_draw < 0.0045)] = "sex_interaction"
    labels[is_x & (x_draw >= 0.0045) & (x_draw < 0.010)] = "female_biased"
    labels[is_x & (x_draw >= 0.010) & (x_draw < 0.014)] = "male_biased"

    is_shared = labels == "shared"
    is_female = labels == "female_biased"
    is_male = labels == "male_biased"
    is_interaction = labels == "sex_interaction"
    is_positive = labels != "null"
    direction = rng.choice([-1, 1], size=size)

    z_total = rng.normal(0, 1, size)
    z_female = rng.normal(0, 1, size)
    z_male = rng.normal(0, 1, size)
    z_interaction = rng.normal(0, 1, size)
    z_total[is_shared] += direction[is_shared] * rng.normal(
        6.4, 1.1, is_shared.sum()
    )
    z_female[is_shared] += direction[is_shared] * rng.normal(
        5.2, 1.0, is_shared.sum()
    )
    z_male[is_shared] += direction[is_shared] * rng.normal(
        5.0, 1.0, is_shared.sum()
    )
    z_female[is_female] += direction[is_female] * rng.normal(
        6.1, 1.0, is_female.sum()
    )
    z_male[is_female] += direction[is_female] * rng.normal(
        1.0, 0.8, is_female.sum()
    )
    z_total[is_female] += direction[is_female] * rng.normal(
        3.5, 1.0, is_female.sum()
    )
    z_interaction[is_female] += direction[is_female] * rng.normal(
        5.0, 0.9, is_female.sum()
    )
    z_male[is_male] += direction[is_male] * rng.normal(5.9, 1.0, is_male.sum())
    z_female[is_male] += direction[is_male] * rng.normal(
        1.0, 0.8, is_male.sum()
    )
    z_total[is_male] += direction[is_male] * rng.normal(
        3.3, 1.0, is_male.sum()
    )
    z_interaction[is_male] -= direction[is_male] * rng.normal(
        4.8, 0.9, is_male.sum()
    )
    z_female[is_interaction] += direction[is_interaction] * rng.normal(
        4.3, 1.0, is_interaction.sum()
    )
    z_male[is_interaction] -= direction[is_interaction] * rng.normal(
        4.0, 1.0, is_interaction.sum()
    )
    z_interaction[is_interaction] += direction[is_interaction] * rng.normal(
        6.5, 1.1, is_interaction.sum()
    )

    se_total = 1 / np.sqrt(150_000 * 2 * maf * (1 - maf))
    se_female = 1 / np.sqrt(85_000 * 2 * maf * (1 - maf))
    se_male = 1 / np.sqrt(65_000 * 2 * maf * (1 - maf))
    se_interaction = np.sqrt(se_female**2 + se_male**2)
    beta_total = z_total * se_total
    beta_female = z_female * se_female
    beta_male = z_male * se_male
    beta_interaction = beta_female - beta_male
    z_interaction = beta_interaction / se_interaction + 0.35 * z_interaction
    p_total = p_from_z(z_total)
    p_female = p_from_z(z_female)
    p_male = p_from_z(z_male)
    p_interaction = p_from_z(z_interaction)

    pip_background = rng.beta(0.3, 30, size)
    pip_total = np.clip(
        pip_background + is_shared * rng.beta(2, 4, size) * 0.75, 0, 1
    )
    pip_female = np.clip(
        pip_background + is_female * rng.beta(2.2, 4, size) * 0.82, 0, 1
    )
    pip_male = np.clip(
        pip_background + is_male * rng.beta(2.2, 4, size) * 0.82, 0, 1
    )
    pip_interaction = np.clip(
        pip_background + is_interaction * rng.beta(2.2, 3.8, size) * 0.88,
        0,
        1,
    )

    tissue, cell, mechanism = expected_targets(labels, is_x, rng)
    is_par = is_x & (
        ((bp >= 10_000) & (bp <= 2_781_479))
        | ((bp >= 155_701_382) & (bp <= 156_030_895))
    )
    is_x_nonpar = is_x & (~is_par)
    xci_draw = rng.random(size)
    xci_escape = is_x_nonpar & (xci_draw < 0.18)
    xci_subject = is_x_nonpar & (xci_draw >= 0.18) & (xci_draw < 0.78)
    xci_variable = is_x_nonpar & (xci_draw >= 0.78)

    gene_pool = np.array(
        [f"GENE{index:05d}" for index in range(1, 50_000)], dtype=object
    )
    causal_gene = rng.choice(gene_pool, size=size)
    causal_gene[is_shared] = rng.choice(KNOWN_AD_GENES, is_shared.sum())
    causal_gene[is_female] = rng.choice(
        FEMALE_GENES + KNOWN_AD_GENES, is_female.sum()
    )
    causal_gene[is_male] = rng.choice(
        MALE_GENES + KNOWN_AD_GENES, is_male.sum()
    )
    causal_gene[is_x] = rng.choice(
        X_GENES + list(gene_pool[:1000]), is_x.sum()
    )

    data: dict[str, object] = {
        "variant_id": variant_id,
        "rsid": rsid,
        "CHR": chromosome,
        "BP": bp,
        "EA": "G",
        "OA": "A",
        "nearest_gene": causal_gene,
        "MAF": maf,
        "INFO": info,
        "LD_score": ld_score,
        "distance_to_tss_kb": distance_to_tss,
        "sex_label": labels,
        "sex_label_id": [CLASS_TO_ID[label] for label in labels],
        "is_positive": is_positive.astype(int),
        "sex_context_id": np.where(is_female, 1, np.where(is_male, 0, 2)),
        "target_gene_label": np.zeros(size, dtype=int),
        "top_tissue_label": tissue,
        "top_cell_label": cell,
        "mechanism_label": mechanism,
        "Z_total": z_total,
        "Z_female": z_female,
        "Z_male": z_male,
        "Z_interaction": z_interaction,
        "P_total": p_total,
        "P_female": p_female,
        "P_male": p_male,
        "P_interaction": p_interaction,
        "neglog10P_total": neglog10(p_total),
        "neglog10P_female": neglog10(p_female),
        "neglog10P_male": neglog10(p_male),
        "neglog10P_interaction": neglog10(p_interaction),
        "SE_total": se_total,
        "SE_female": se_female,
        "SE_male": se_male,
        "SE_interaction": se_interaction,
        "Beta_total": beta_total,
        "Beta_female": beta_female,
        "Beta_male": beta_male,
        "Beta_interaction": beta_interaction,
        "sex_delta_beta": beta_female - beta_male,
        "abs_sex_delta_beta": np.abs(beta_female - beta_male),
        "PIP_total": pip_total,
        "PIP_female": pip_female,
        "PIP_male": pip_male,
        "PIP_interaction": pip_interaction,
        "deepsea_delta": np.clip(
            rng.normal(0, 0.5, size)
            + is_positive * rng.normal(1, 0.4, size),
            -3,
            5,
        ),
        "enformer_brain_delta": np.clip(
            rng.normal(0, 0.5, size)
            + (is_shared | is_female | is_interaction)
            * rng.normal(1.2, 0.4, size),
            -3,
            5,
        ),
        "sei_regulatory_score": np.clip(
            rng.normal(0, 0.5, size)
            + (is_shared | is_male | is_female) * rng.normal(0.9, 0.4, size),
            -3,
            5,
        ),
        "cadd_score": np.clip(
            rng.normal(10, 4, size)
            + is_positive * rng.normal(5, 2, size),
            0,
            50,
        ),
        "conservation_score": np.clip(
            rng.beta(1.2, 3, size) + is_positive * 0.15, 0, 1
        ),
        "chr_type": np.where(
            is_x, np.where(is_par, "PAR", "X_nonPAR"), "autosome"
        ),
        "xci_status": np.where(
            is_par,
            "PAR",
            np.where(
                xci_escape,
                "escape",
                np.where(
                    xci_subject,
                    "subject",
                    np.where(xci_variable, "variable", "autosome"),
                ),
            ),
        ),
        "x_is_chrX": is_x.astype(int),
        "x_is_nonPAR": is_x_nonpar.astype(int),
        "xci_escape": xci_escape.astype(int),
        "xci_subject": xci_subject.astype(int),
        "xci_variable": xci_variable.astype(int),
    }

    for candidate in range(6):
        boost = is_positive.astype(float) if candidate == 0 else np.zeros(size)
        names = (
            causal_gene.copy()
            if candidate == 0
            else rng.choice(gene_pool, size=size).astype(object)
        )
        prefix = f"gene{candidate}_"
        data[prefix + "name"] = names
        data[prefix + "distance_kb"] = np.clip(
            np.exp(rng.normal(4, 1, size)) - boost * rng.normal(35, 10, size),
            0,
            10_000,
        )
        data[prefix + "eqtl_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            2, 0.7, size
        )
        data[prefix + "sqtl_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            1.5, 0.6, size
        )
        data[prefix + "pqtl_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            1.2, 0.6, size
        )
        data[prefix + "twas_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            2.2, 0.7, size
        )
        data[prefix + "coloc_pph4"] = np.clip(
            rng.beta(0.8, 10, size) + boost * rng.beta(2.2, 3.2, size) * 0.75,
            0,
            1,
        )
        data[prefix + "brain_expression"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(1.2, 0.5, size),
            -4,
            8,
        )
        data[prefix + "ad_de_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            1.8, 0.7, size
        )
        data[prefix + "constraint_score"] = np.clip(
            rng.beta(1.5, 4, size) + boost * 0.1, 0, 1
        )

    for tissue_index in range(len(TISSUES)):
        boost = (tissue == tissue_index).astype(float)
        prefix = f"tissue{tissue_index}_"
        data[prefix + "eqtl_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            2, 0.7, size
        )
        data[prefix + "sqtl_z"] = rng.normal(0, 1, size) + boost * rng.normal(
            1.3, 0.6, size
        )
        data[prefix + "coloc_pph4"] = np.clip(
            rng.beta(0.8, 10, size) + boost * rng.beta(2, 3, size) * 0.75,
            0,
            1,
        )
        data[prefix + "enhancer"] = np.clip(
            rng.beta(1, 6, size) + boost * rng.beta(2.2, 4, size) * 0.7,
            0,
            1,
        )
        data[prefix + "ldsc_enrichment"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(1.4, 0.5, size),
            -4,
            6,
        )
        data[prefix + "expression"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(1.3, 0.6, size),
            -4,
            8,
        )

    for cell_index in range(len(CELL_TYPES)):
        boost = (cell == cell_index).astype(float)
        prefix = f"cell{cell_index}_"
        data[prefix + "marker_score"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(2, 0.6, size), -4, 8
        )
        data[prefix + "ad_de_score"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(1.7, 0.6, size),
            -4,
            8,
        )
        data[prefix + "enhancer_score"] = np.clip(
            rng.beta(1, 6, size) + boost * rng.beta(2.2, 3.5, size) * 0.75,
            0,
            1,
        )
        data[prefix + "expression_specificity"] = np.clip(
            rng.normal(0, 1, size) + boost * rng.normal(1.4, 0.5, size),
            -4,
            8,
        )

    data["coloc_brain_pph4"] = np.maximum.reduce(
        [
            data["tissue0_coloc_pph4"],
            data["tissue1_coloc_pph4"],
            data["tissue2_coloc_pph4"],
        ]
    )
    data["coloc_blood_pph4"] = data["tissue3_coloc_pph4"]
    data["coloc_immune_pph4"] = np.maximum(
        data["tissue4_coloc_pph4"], data["tissue5_coloc_pph4"]
    )
    return pd.DataFrame(data)


def main() -> int:
    args = parse_args()
    if args.n_variants <= 0 or args.chunk_size <= 0:
        raise ValueError("--n-variants and --chunk-size must be positive")

    outdir = Path(args.outdir)
    ensure_dir(outdir)
    partition_dir = outdir / "simulated_by_chr"
    ensure_dir(partition_dir)
    for old_part in partition_dir.glob("chr=*/part_*.parquet"):
        old_part.unlink()

    output_path = outdir / "simulated_variants.parquet"
    output_path.unlink(missing_ok=True)
    rng = np.random.default_rng(args.seed)
    label_counts = Counter()
    writer = None
    n_chunks = int(np.ceil(args.n_variants / args.chunk_size))

    try:
        for index in tqdm(range(n_chunks), desc="simulate-v3"):
            start = index * args.chunk_size
            size = min(args.chunk_size, args.n_variants - start)
            frame = simulate_chunk(start, size, rng)
            label_counts.update(frame["sex_label"].value_counts().to_dict())

            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema)
            writer.write_table(table)

            for chromosome, subset in frame.groupby("CHR", sort=False):
                chromosome_dir = partition_dir / f"chr={chromosome}"
                ensure_dir(chromosome_dir)
                subset.to_parquet(
                    chromosome_dir / f"part_{index:05d}.parquet", index=False
                )
    finally:
        if writer is not None:
            writer.close()

    label_table = pd.DataFrame(
        sorted(label_counts.items()), columns=["sex_label", "n"]
    )
    label_table.to_csv(outdir / "label_distribution.csv", index=False)
    print(f"Wrote {args.n_variants:,} variants to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
