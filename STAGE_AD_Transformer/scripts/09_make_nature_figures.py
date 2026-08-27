#!/usr/bin/env python3
"""Step 09: generate high-resolution summary figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


CHROMOSOMES = [str(index) for index in range(1, 23)] + ["X"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate SexReg-AD figures.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--ablation", required=True)
    parser.add_argument("--outdir", default="results/figures")
    parser.add_argument("--max-points", type=int, default=1_000_000)
    return parser.parse_args()


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=450, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def manhattan(frame: pd.DataFrame, p_column: str, title: str, path: Path) -> None:
    data = frame[["CHR", "BP", p_column]].dropna().copy()
    data = data[data["CHR"].astype(str).isin(CHROMOSOMES)]
    if data.empty:
        print(f"Skip {path.name}: no valid chromosome data")
        return
    data["CHR"] = pd.Categorical(
        data["CHR"].astype(str), categories=CHROMOSOMES, ordered=True
    )
    data = data.sort_values(["CHR", "BP"])

    x_segments = []
    ticks = []
    labels = []
    offset = 0
    for chromosome in CHROMOSOMES:
        subset = data[data["CHR"] == chromosome]
        if subset.empty:
            continue
        positions = subset["BP"].to_numpy() + offset
        x_segments.append(pd.Series(positions, index=subset.index))
        ticks.append(positions.mean())
        labels.append(chromosome)
        offset = positions.max() + 5_000_000

    data["x"] = pd.concat(x_segments).sort_index()
    data["y"] = -np.log10(np.clip(data[p_column], 1e-300, 1))
    plt.figure(figsize=(7.4, 2.5))
    codes = data["CHR"].cat.codes.to_numpy()
    for parity in (0, 1):
        mask = codes % 2 == parity
        plt.scatter(
            data.loc[mask, "x"],
            data.loc[mask, "y"],
            s=1.4,
            alpha=0.65,
            linewidths=0,
        )
    plt.axhline(-np.log10(5e-8), linestyle="--", linewidth=0.8)
    plt.axhline(-np.log10(1e-6), linestyle=":", linewidth=0.8)
    plt.xticks(ticks, labels)
    plt.xlabel("Chromosome")
    plt.ylabel(r"$-log_{10}(P)$")
    plt.title(title)
    save_figure(path)


def plot_benchmark(path: str, outdir: Path) -> None:
    metrics = pd.read_csv(path).sort_values("AUPRC", ascending=True)
    plt.figure(figsize=(5.5, 3.4))
    plt.barh(metrics["model"], metrics["AUPRC"])
    plt.xlabel("AUPRC")
    plt.title("Benchmark against traditional methods")
    save_figure(outdir / "fig5_model_benchmark.png")


def plot_ablation(path: str, outdir: Path) -> None:
    metrics = pd.read_csv(path)
    full_rows = metrics[metrics["ablation"] == "full"]
    if full_rows.empty:
        raise ValueError("Ablation metrics do not contain a 'full' row")
    full_auprc = full_rows["AUPRC"].iloc[0]
    metrics = metrics[metrics["ablation"] != "full"].copy()
    metrics["drop"] = full_auprc - metrics["AUPRC"]
    metrics = metrics.sort_values("drop")
    plt.figure(figsize=(5.3, 3.4))
    plt.barh(metrics["ablation"], metrics["drop"])
    plt.xlabel("AUPRC decrease vs full model")
    plt.title("Inference-time ablation analysis")
    save_figure(outdir / "fig6_ablation.png")


def plot_prediction_summaries(predictions: pd.DataFrame, outdir: Path) -> None:
    sample = (
        predictions.sample(200_000, random_state=1)
        if len(predictions) > 200_000
        else predictions
    )
    plt.figure(figsize=(3.4, 3.2))
    plt.scatter(
        sample["female_context_score"],
        sample["male_context_score"],
        s=3,
        alpha=0.3,
        linewidths=0,
    )
    plt.xlabel("Female-context score")
    plt.ylabel("Male-context score")
    plt.title("Sex-context scores")
    save_figure(outdir / "fig7_context_score_scatter.png")

    table = pd.crosstab(
        predictions["predicted_class"],
        predictions["predicted_cell_type"],
        normalize="index",
    )
    plt.figure(figsize=(6, 2.8))
    image = plt.imshow(table.to_numpy(), aspect="auto")
    plt.colorbar(image, fraction=0.025, pad=0.02, label="Row fraction")
    plt.yticks(range(len(table.index)), table.index)
    plt.xticks(range(len(table.columns)), table.columns, rotation=45, ha="right")
    plt.title("Predicted cell-type distribution by sex class")
    save_figure(outdir / "fig8_celltype_heatmap.png")


def plot_chromosome_x(data: pd.DataFrame, outdir: Path) -> None:
    chromosome_x = (
        data[data["CHR"].astype(str) == "X"]
        .sort_values("P_interaction")
        .head(100)
        .copy()
    )
    if chromosome_x.empty:
        print("Skip fig9_chrX_lollipop: no chromosome X rows")
        return
    chromosome_x["y"] = -np.log10(
        np.clip(chromosome_x["P_interaction"], 1e-300, 1)
    )
    plt.figure(figsize=(7.2, 2.6))
    plt.vlines(chromosome_x["BP"], 0, chromosome_x["y"], linewidth=0.7, alpha=0.6)
    plt.scatter(chromosome_x["BP"], chromosome_x["y"], s=12, alpha=0.85)
    plt.axhline(-np.log10(5e-8), linestyle="--", linewidth=0.8)
    plt.xlabel("Position on chromosome X")
    plt.ylabel(r"$-log_{10}(P_{interaction})$")
    plt.title("Chromosome X SNP × sex interaction signals")
    save_figure(outdir / "fig9_chrX_lollipop.png")


def plot_network(predictions: pd.DataFrame, outdir: Path) -> None:
    top = predictions.sort_values("sexreg_ad_score", ascending=False).head(25)
    graph = nx.Graph()
    for _, row in top.iterrows():
        chain = [
            row["rsid"],
            row["predicted_target_gene"],
            row["predicted_tissue"],
            row["predicted_cell_type"],
            row["predicted_mechanism"],
        ]
        graph.add_nodes_from(chain)
        graph.add_edges_from(zip(chain[:-1], chain[1:]))
    positions = nx.spring_layout(graph, seed=4, k=0.85)
    plt.figure(figsize=(7, 5.2))
    nx.draw_networkx_edges(graph, positions, alpha=0.25, width=0.7)
    nx.draw_networkx_nodes(graph, positions, node_size=90, alpha=0.9)
    nx.draw_networkx_labels(graph, positions, font_size=5)
    plt.axis("off")
    plt.title("Top-ranked variant-gene-tissue-cell-mechanism network")
    save_figure(outdir / "fig10_network.png")


def main() -> int:
    args = parse_args()
    if args.max_points <= 0:
        raise ValueError("--max-points must be positive")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    configure_style()

    data = pd.read_parquet(args.data)
    figure_data = (
        data.sample(args.max_points, random_state=1)
        if len(data) > args.max_points
        else data
    )
    manhattan(
        figure_data,
        "P_total",
        "Sex-combined AD GWAS",
        outdir / "fig1_manhattan_total.png",
    )
    manhattan(
        figure_data,
        "P_female",
        "Female-stratified AD GWAS",
        outdir / "fig2_manhattan_female.png",
    )
    manhattan(
        figure_data,
        "P_male",
        "Male-stratified AD GWAS",
        outdir / "fig3_manhattan_male.png",
    )
    manhattan(
        figure_data,
        "P_interaction",
        "SNP × sex interaction GWAS",
        outdir / "fig4_manhattan_interaction.png",
    )
    plot_benchmark(args.benchmark, outdir)
    plot_ablation(args.ablation, outdir)

    predictions = pd.read_csv(args.predictions)
    plot_prediction_summaries(predictions, outdir)
    plot_chromosome_x(data, outdir)
    plot_network(predictions, outdir)
    print(f"Figures written to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
