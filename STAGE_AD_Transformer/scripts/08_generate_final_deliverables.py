#!/usr/bin/env python3
"""Step 08: assemble publication-oriented result tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.constants import KNOWN_AD_GENES
from sexreg_ad.utils import ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assemble final result tables.")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--baseline-metrics", required=True)
    parser.add_argument("--transformer-metrics", required=True)
    parser.add_argument("--ablation-metrics", required=True)
    parser.add_argument(
        "--model", default="results/checkpoints/sexreg_ad_finetuned.pt"
    )
    parser.add_argument("--outdir", default="results/final_outputs")
    parser.add_argument("--top-n", type=int, default=500)
    parser.add_argument("--locus-window-bp", type=int, default=1_000_000)
    return parser.parse_args()


def save_selected(
    frame: pd.DataFrame,
    columns: list[str],
    path: Path,
    *,
    sort_by: str | None = None,
) -> None:
    selected = frame[[column for column in columns if column in frame]].copy()
    if sort_by is not None:
        selected = selected.sort_values(sort_by)
    selected.to_csv(path, index=False)


def main() -> int:
    args = parse_args()
    if args.top_n <= 0 or args.locus_window_bp <= 0:
        raise ValueError("--top-n and --locus-window-bp must be positive")
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    predictions = pd.read_csv(args.predictions)
    predictions["locus_id"] = (
        "chr"
        + predictions["CHR"].astype(str)
        + ":bin"
        + (predictions["BP"].astype(float) // args.locus_window_bp)
        .astype(int)
        .astype(str)
    )
    predictions["ad_risk_priority_rank"] = predictions[
        "sexreg_ad_score"
    ].rank(ascending=False, method="min").astype(int)
    predictions["priority_tier"] = pd.cut(
        predictions["sexreg_ad_score"].rank(pct=True),
        bins=[0, 0.9, 0.98, 0.995, 1],
        labels=["Tier 4", "Tier 3", "Tier 2", "Tier 1"],
        include_lowest=True,
    )

    pd.DataFrame(
        [{"deliverable": "SexReg-AD Transformer v3 model", "path": args.model}]
    ).to_csv(outdir / "00_model_registry.csv", index=False)

    priority_columns = [
        "variant_id",
        "rsid",
        "CHR",
        "BP",
        "locus_id",
        "nearest_gene",
        "predicted_target_gene",
        "sexreg_ad_score",
        "ad_risk_priority_rank",
        "priority_tier",
        "predicted_class",
        "female_context_score",
        "male_context_score",
        "delta_sex_score",
        "predicted_tissue",
        "predicted_cell_type",
        "predicted_mechanism",
        "P_total",
        "P_female",
        "P_male",
        "P_interaction",
        "chr_type",
        "xci_status",
    ]
    save_selected(
        predictions,
        priority_columns,
        outdir / "01_variant_ad_risk_priority_table.csv",
        sort_by="ad_risk_priority_rank",
    )

    sex_columns = [
        "variant_id",
        "rsid",
        "CHR",
        "BP",
        "locus_id",
        "nearest_gene",
        "predicted_class",
        "sexreg_ad_score",
        "female_context_score",
        "male_context_score",
        "delta_sex_score",
        *[column for column in predictions if column.startswith("prob_class_")],
    ]
    save_selected(
        predictions,
        sex_columns,
        outdir / "02_sex_class_prediction_table.csv",
    )

    lead = (
        predictions.sort_values("sexreg_ad_score", ascending=False)
        .groupby("locus_id", as_index=False)
        .head(1)
        .copy()
    )
    save_selected(
        lead,
        [
            "locus_id",
            "CHR",
            "BP",
            "rsid",
            "nearest_gene",
            "predicted_target_gene",
            "sexreg_ad_score",
            "predicted_class",
            "P_total",
            "P_female",
            "P_male",
            "P_interaction",
        ],
        outdir / "03_top_target_genes_per_locus.csv",
    )
    save_selected(
        lead,
        ["locus_id", "rsid", "predicted_class", "predicted_tissue", "sexreg_ad_score"],
        outdir / "04_top_tissues_per_locus.csv",
    )
    save_selected(
        lead,
        [
            "locus_id",
            "rsid",
            "predicted_class",
            "predicted_cell_type",
            "sexreg_ad_score",
        ],
        outdir / "05_top_cell_types_per_locus.csv",
    )

    ablations = pd.read_csv(args.ablation_metrics)
    full_rows = ablations[ablations["ablation"] == "full"]
    if full_rows.empty:
        raise ValueError("Ablation metrics do not contain a 'full' row")
    full = full_rows.iloc[0]
    attribution_rows = []
    for _, row in ablations.iterrows():
        if row["ablation"] != "full":
            attribution_rows.append(
                {
                    "feature_domain": row["ablation"],
                    "delta_AUPRC_vs_full": full["AUPRC"] - row["AUPRC"],
                    "delta_AUROC_vs_full": full["AUROC"] - row["AUROC"],
                }
            )
    pd.DataFrame(attribution_rows).sort_values(
        "delta_AUPRC_vs_full", ascending=False
    ).to_csv(outdir / "06_model_attribution_by_domain.csv", index=False)

    benchmark = pd.concat(
        [
            pd.read_csv(args.baseline_metrics),
            pd.read_csv(args.transformer_metrics),
        ],
        ignore_index=True,
    ).sort_values("AUPRC", ascending=False)
    benchmark.to_csv(
        outdir / "07_benchmark_against_traditional_methods.csv", index=False
    )

    known_genes = set(KNOWN_AD_GENES)
    lead["is_known_ad_gene"] = lead["nearest_gene"].isin(known_genes) | lead[
        "predicted_target_gene"
    ].isin(known_genes)
    novel = (
        lead[~lead["is_known_ad_gene"]]
        .sort_values("sexreg_ad_score", ascending=False)
        .head(args.top_n)
    )
    candidate_columns = [
        "locus_id",
        "rsid",
        "CHR",
        "BP",
        "nearest_gene",
        "predicted_target_gene",
        "sexreg_ad_score",
        "predicted_class",
        "predicted_tissue",
        "predicted_cell_type",
        "predicted_mechanism",
        "P_total",
        "P_female",
        "P_male",
        "P_interaction",
        "chr_type",
        "xci_status",
    ]
    save_selected(
        novel,
        candidate_columns,
        outdir / "08_top_ranked_novel_loci.csv",
    )

    aim3 = lead.sort_values("sexreg_ad_score", ascending=False).head(args.top_n).copy()
    aim3["recommended_validation"] = (
        "sex-stratified replication + SNP×sex interaction + fine-mapping + "
        "colocalization + TWAS/sTWAS/PWAS + biomarker validation"
    )
    save_selected(
        aim3,
        [
            column
            for column in candidate_columns
            if column not in {"nearest_gene"}
        ]
        + ["recommended_validation"],
        outdir / "09_aim3_candidate_mechanisms.csv",
    )

    manifest = pd.DataFrame(
        [
            [1, "SexReg-AD Transformer model", "00_model_registry.csv"],
            [2, "Variant risk priority", "01_variant_ad_risk_priority_table.csv"],
            [3, "Sex-class prediction", "02_sex_class_prediction_table.csv"],
            [4, "Top target genes", "03_top_target_genes_per_locus.csv"],
            [5, "Top tissues", "04_top_tissues_per_locus.csv"],
            [6, "Top cell types", "05_top_cell_types_per_locus.csv"],
            [7, "Explainability", "06_model_attribution_by_domain.csv"],
            [8, "Benchmark", "07_benchmark_against_traditional_methods.csv"],
            [9, "Novel loci", "08_top_ranked_novel_loci.csv"],
            [10, "Aim 3 mechanisms", "09_aim3_candidate_mechanisms.csv"],
        ],
        columns=["item", "deliverable", "file"],
    )
    manifest.to_csv(outdir / "MANIFEST_final_deliverables.csv", index=False)
    print(manifest.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
