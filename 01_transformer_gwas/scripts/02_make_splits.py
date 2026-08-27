#!/usr/bin/env python3
"""Step 02: create chromosome-held-out train, validation, and test splits."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as pads
import pyarrow.parquet as pq

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.utils import ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream a Parquet dataset into chromosome-aware splits."
    )
    parser.add_argument("--input", required=True, help="Input Parquet file/directory.")
    parser.add_argument("--outdir", default="data/splits")
    parser.add_argument(
        "--holdout-chrs",
        nargs="+",
        default=["8", "16", "X"],
        help="Chromosomes assigned exclusively to the test split.",
    )
    parser.add_argument("--valid-frac", type=float, default=0.15)
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260618)
    return parser.parse_args()


def iter_batches(path: Path, batch_size: int) -> Iterable[pa.RecordBatch]:
    if path.is_dir():
        yield from pads.dataset(path, format="parquet").to_batches(batch_size=batch_size)
    else:
        yield from pq.ParquetFile(path).iter_batches(batch_size=batch_size)


def main() -> int:
    args = parse_args()
    if not 0 < args.valid_frac < 1:
        raise ValueError("--valid-frac must be between 0 and 1")
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    output_paths = {
        name: outdir / f"{name}.parquet" for name in ("train", "valid", "test")
    }
    for path in output_paths.values():
        path.unlink(missing_ok=True)

    holdout = {str(chromosome) for chromosome in args.holdout_chrs}
    rng = np.random.default_rng(args.seed)
    writers: dict[str, pq.ParquetWriter] = {}
    stats = {
        name: {"n": 0, "positives": 0, "chromosomes": set()}
        for name in output_paths
    }

    try:
        for batch in iter_batches(input_path, args.batch_size):
            frame = batch.to_pandas()
            if "CHR" not in frame or "is_positive" not in frame:
                raise ValueError("Input must contain CHR and is_positive columns")
            is_test = frame["CHR"].astype(str).isin(holdout).to_numpy()
            is_valid = (~is_test) & (rng.random(len(frame)) < args.valid_frac)
            masks = {
                "train": (~is_test) & (~is_valid),
                "valid": is_valid,
                "test": is_test,
            }

            for name, mask in masks.items():
                subset = frame.loc[mask]
                if subset.empty:
                    continue
                table = pa.Table.from_pandas(subset, preserve_index=False)
                if name not in writers:
                    writers[name] = pq.ParquetWriter(output_paths[name], table.schema)
                writers[name].write_table(table)
                stats[name]["n"] += len(subset)
                stats[name]["positives"] += int(subset["is_positive"].sum())
                stats[name]["chromosomes"].update(
                    subset["CHR"].astype(str).unique().tolist()
                )
    finally:
        for writer in writers.values():
            writer.close()

    empty = [name for name, values in stats.items() if values["n"] == 0]
    if empty:
        raise ValueError(f"Empty split(s): {', '.join(empty)}")

    summary_rows = []
    for name in ("train", "valid", "test"):
        values = stats[name]
        summary_rows.append(
            {
                "split": name,
                "n": values["n"],
                "positive_rate": values["positives"] / values["n"],
                "chromosomes": ",".join(
                    sorted(values["chromosomes"], key=lambda x: (x == "X", int(x) if x.isdigit() else 99))
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "split_summary.csv", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
