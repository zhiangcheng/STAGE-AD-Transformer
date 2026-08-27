# Contributing

Thank you for contributing to STAGE-AD Transformer.

1. Open an issue describing the scientific or software problem.
2. Create a focused branch and keep generated data/checkpoints out of Git.
3. Add or update tests for behavior changes.
4. Run `python -m compileall -q stage_ad scripts run_pipeline.py` and
   `python -m pytest -q`.
5. Explain changes to data assumptions, model outputs, or evaluation protocols
   in the pull request.

Do not commit participant-level data, credentials, access tokens, controlled
datasets, or identifiable metadata. Synthetic examples must be clearly labeled
as synthetic. Scientific claims should state their evidence and limitations.
