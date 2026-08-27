#!/usr/bin/env python3
"""Step 00: report environment readiness."""

from __future__ import annotations

import importlib
import platform
import sys


CORE = [
    ("numpy", "NumPy"),
    ("pandas", "pandas"),
    ("pyarrow", "PyArrow"),
    ("scipy", "SciPy"),
    ("sklearn", "scikit-learn"),
    ("matplotlib", "Matplotlib"),
    ("torch", "PyTorch"),
    ("tqdm", "tqdm"),
    ("networkx", "NetworkX"),
]
OPTIONAL = [("xgboost", "XGBoost"), ("lightgbm", "LightGBM")]


def import_status(module_name: str, display_name: str, *, required: bool) -> bool:
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, "__version__", "version unknown")
        label = "OK" if required else "OK optional"
        print(f"[{label}] {display_name}: {version}")
        return True
    except Exception as exc:
        label = "MISSING" if required else "MISSING optional"
        print(f"[{label}] {display_name}: {exc}")
        return not required


def main() -> int:
    print(f"Python: {platform.python_version()} ({sys.executable})")
    if sys.version_info < (3, 10):
        print("[ERROR] Python 3.10 or newer is required.")
        return 1

    statuses = [
        import_status(module, name, required=True) for module, name in CORE
    ]
    ready = all(statuses)
    for module, name in OPTIONAL:
        import_status(module, name, required=False)

    if ready:
        import torch

        device = "CUDA" if torch.cuda.is_available() else "CPU"
        print(f"Compute device detected: {device}")
        print("Dependency check completed.")
        return 0

    print("Install the project first: python -m pip install -e .")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
