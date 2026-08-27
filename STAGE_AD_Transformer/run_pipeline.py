#!/usr/bin/env python3
"""Run the SexReg-AD numbered workflow with a reproducible profile."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROFILE_PATH = ROOT / "config" / "pipeline_profiles.json"

STEP_NAMES = {
    0: "check dependencies",
    1: "prepare synthetic data",
    2: "create chromosome-aware splits",
    3: "pretrain regulatory representations",
    4: "fine-tune the multi-task model",
    5: "train and evaluate baselines",
    6: "generate predictions",
    7: "run inference-time ablations",
    8: "assemble final result tables",
    9: "generate figures",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Execute the numbered SexReg-AD pipeline in order."
    )
    parser.add_argument(
        "--profile",
        default="smoke",
        help="Profile name from config/pipeline_profiles.json (default: smoke).",
    )
    parser.add_argument(
        "--from-step",
        type=int,
        choices=range(10),
        default=0,
        help="First step to run, inclusive (default: 0).",
    )
    parser.add_argument(
        "--to-step",
        type=int,
        choices=range(10),
        default=9,
        help="Last step to run, inclusive (default: 9).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    parser.add_argument(
        "--list-profiles",
        action="store_true",
        help="List available profiles and exit.",
    )
    return parser.parse_args()


def load_profiles() -> dict[str, dict[str, object]]:
    with PROFILE_PATH.open(encoding="utf-8") as handle:
        profiles = json.load(handle)
    if not isinstance(profiles, dict) or not profiles:
        raise ValueError(f"No profiles found in {PROFILE_PATH}")
    return profiles


def build_commands(
    config: dict[str, object], python: str
) -> dict[int, list[str]]:
    scripts = ROOT / "scripts"
    processed = ROOT / "data" / "processed" / "simulated_variants.parquet"
    splits = ROOT / "data" / "splits"
    results = ROOT / "results"
    metrics = results / "metrics"
    checkpoint_dir = results / "checkpoints"
    predictions = results / "predictions" / "predictions.csv"
    final_outputs = results / "final_outputs"

    common_model_args = [
        "--batch-size",
        str(config["batch_size"]),
        "--max-train",
        str(config["max_train"]),
        "--d-model",
        str(config["d_model"]),
        "--n-heads",
        str(config["n_heads"]),
        "--n-layers",
        str(config["n_layers"]),
        "--seed",
        str(config["seed"]),
    ]

    return {
        0: [python, str(scripts / "00_check_dependencies.py")],
        1: [
            python,
            str(scripts / "01_simulate_data.py"),
            "--n-variants",
            str(config["n_variants"]),
            "--chunk-size",
            str(config["chunk_size"]),
            "--outdir",
            str(processed.parent),
            "--seed",
            str(config["seed"]),
        ],
        2: [
            python,
            str(scripts / "02_make_splits.py"),
            "--input",
            str(processed),
            "--outdir",
            str(splits),
            "--holdout-chrs",
            *[str(value) for value in config["holdout_chrs"]],
            "--valid-frac",
            str(config["valid_frac"]),
            "--seed",
            str(config["seed"]),
        ],
        3: [
            python,
            str(scripts / "03_pretrain_regulatory.py"),
            "--train",
            str(splits / "train.parquet"),
            "--valid",
            str(splits / "valid.parquet"),
            "--test",
            str(splits / "test.parquet"),
            "--outdir",
            str(results),
            "--epochs",
            str(config["pretrain_epochs"]),
            *common_model_args,
        ],
        4: [
            python,
            str(scripts / "04_finetune_multitask.py"),
            "--train",
            str(splits / "train.parquet"),
            "--valid",
            str(splits / "valid.parquet"),
            "--test",
            str(splits / "test.parquet"),
            "--pretrained",
            str(checkpoint_dir / "sexreg_ad_pretrained.pt"),
            "--outdir",
            str(results),
            "--epochs",
            str(config["finetune_epochs"]),
            *common_model_args,
        ],
        5: [
            python,
            str(scripts / "05_train_baselines.py"),
            "--train",
            str(splits / "train.parquet"),
            "--test",
            str(splits / "test.parquet"),
            "--outdir",
            str(metrics),
            "--max-train",
            str(config["max_train"]),
            "--seed",
            str(config["seed"]),
        ],
        6: [
            python,
            str(scripts / "06_predict.py"),
            "--model",
            str(checkpoint_dir / "sexreg_ad_finetuned.pt"),
            "--input",
            str(splits / "test.parquet"),
            "--outdir",
            str(results),
            "--batch-size",
            str(config["predict_batch_size"]),
        ],
        7: [
            python,
            str(scripts / "07_run_ablations.py"),
            "--model",
            str(checkpoint_dir / "sexreg_ad_finetuned.pt"),
            "--test",
            str(splits / "test.parquet"),
            "--outdir",
            str(metrics),
            "--batch-size",
            str(config["predict_batch_size"]),
        ],
        8: [
            python,
            str(scripts / "08_generate_final_deliverables.py"),
            "--predictions",
            str(predictions),
            "--baseline-metrics",
            str(metrics / "baseline_metrics.csv"),
            "--transformer-metrics",
            str(metrics / "transformer_test_metrics.csv"),
            "--ablation-metrics",
            str(metrics / "ablation_metrics.csv"),
            "--model",
            str(checkpoint_dir / "sexreg_ad_finetuned.pt"),
            "--outdir",
            str(final_outputs),
            "--top-n",
            str(config["top_n"]),
        ],
        9: [
            python,
            str(scripts / "09_make_nature_figures.py"),
            "--data",
            str(processed),
            "--predictions",
            str(predictions),
            "--benchmark",
            str(final_outputs / "07_benchmark_against_traditional_methods.csv"),
            "--ablation",
            str(metrics / "ablation_metrics.csv"),
            "--outdir",
            str(results / "figures"),
            "--max-points",
            str(config["max_figure_points"]),
        ],
    }


def display_command(command: Sequence[str]) -> str:
    return shlex.join(str(part) for part in command)


def main() -> int:
    args = parse_args()
    profiles = load_profiles()

    if args.list_profiles:
        print("Available profiles:")
        for name in profiles:
            print(f"  - {name}")
        return 0

    if args.profile not in profiles:
        choices = ", ".join(profiles)
        raise SystemExit(f"Unknown profile {args.profile!r}. Choose one of: {choices}")
    if args.from_step > args.to_step:
        raise SystemExit("--from-step must be less than or equal to --to-step")

    commands = build_commands(profiles[args.profile], sys.executable)
    print(f"Profile: {args.profile}")
    print(f"Steps: {args.from_step}-{args.to_step}")

    for step in range(args.from_step, args.to_step + 1):
        command = commands[step]
        print(f"\n[{step:02d}/09] {STEP_NAMES[step]}", flush=True)
        print(f"$ {display_command(command)}", flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)

    print("\nPipeline completed." if not args.dry_run else "\nDry run completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
