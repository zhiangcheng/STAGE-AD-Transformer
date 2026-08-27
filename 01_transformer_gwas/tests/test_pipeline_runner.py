import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pipeline_dry_run_lists_numbered_steps():
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "run_pipeline.py"),
            "--profile",
            "smoke",
            "--from-step",
            "0",
            "--to-step",
            "2",
            "--dry-run",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "[00/09] check dependencies" in result.stdout
    assert "[01/09] prepare synthetic data" in result.stdout
    assert "[02/09] create chromosome-aware splits" in result.stdout
    assert "Dry run completed." in result.stdout
