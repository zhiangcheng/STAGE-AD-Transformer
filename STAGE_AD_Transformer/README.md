# SexReg-AD Transformer

Workflow | [Data schema](docs/DATA_SCHEMA.md) | [Output reference](docs/OUTPUTS.md)

SexReg-AD Transformer is a research pipeline for sex-aware prioritization of
Alzheimer's disease (AD) variants, candidate genes, tissues, cell types, and
regulatory mechanisms. It combines sex-stratified GWAS features with QTL,
fine-mapping, sequence, tissue, and single-cell evidence in a hierarchical
multi-task Transformer.

> **Research-use warning:** the bundled generator creates synthetic data and
> synthetic labels for software validation only. Outputs from the demo are not
> biological discoveries and must replace your own harmonized real data.

## Workflow

| Step | Script | Main input | Main output |
|---:|---|---|---|
| 00 | `scripts/00_check_dependencies.py` | Python environment | dependency report |
| 01 | `scripts/01_simulate_data.py` | profile/CLI arguments | `data/processed/simulated_variants.parquet` |
| 02 | `scripts/02_make_splits.py` | processed Parquet | train/valid/test Parquet files |
| 03 | `scripts/03_pretrain_regulatory.py` | three splits | pretrained checkpoint |
| 04 | `scripts/04_finetune_multitask.py` | splits + pretrained checkpoint | fine-tuned checkpoint + test metrics |
| 05 | `scripts/05_train_baselines.py` | train/test splits | baseline metrics |
| 06 | `scripts/06_predict.py` | fine-tuned checkpoint + input data | prediction CSV/Parquet |
| 07 | `scripts/07_run_ablations.py` | checkpoint + test split | ablation metrics |
| 08 | `scripts/08_generate_final_deliverables.py` | predictions + metrics | result tables |
| 09 | `scripts/09_make_nature_figures.py` | data + results | PNG/PDF figures |

## Installation

Python 3.10 or 3.11 is recommended. GPU support is optional; PyTorch will use
CUDA when it is available.

### Conda

```bash
conda env create -f environment.yml
conda activate sexreg-ad
python -m pip install -e .
```

### venv + pip

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

Optional XGBoost and LightGBM baselines:

```bash
python -m pip install -e ".[baselines]"
```

Verify the environment:

```bash
python scripts/00_check_dependencies.py
```

## Quick start

The cross-platform runner executes the numbered scripts in order. Start with
the small smoke profile to verify the installation:

```bash
python run_pipeline.py --profile smoke
```

Run the standard synthetic demonstration:

```bash
python run_pipeline.py --profile demo
```

Preview commands without running them:

```bash
python run_pipeline.py --profile demo --dry-run
```

Resume from a selected range after existing outputs have been checked:

```bash
python run_pipeline.py --profile demo --from-step 6 --to-step 9
```

The `large-10m` profile generates GB-scale data and is intended for a server:

```bash
python run_pipeline.py --profile large-10m
```

All profile parameters are visible in
[`config/pipeline_profiles.json`](config/pipeline_profiles.json). CLI commands
printed by the runner can also be copied and executed one at a time.

## Repository layout

```text
.
|-- config/                 Reproducible run profiles and project metadata
|-- docs/                   Step-by-step, schema, and output documentation
|-- scripts/                Numbered pipeline entry points (00-09)
|-- sexreg_ad/              Reusable dataset, model, loss, and metric code
|-- tests/                  Fast unit/smoke tests
|-- run_pipeline.py         Cross-platform workflow runner
|-- environment.yml         Conda environment
|-- pyproject.toml          Installable Python package metadata
`-- requirements*.txt       pip dependency lists
```

Generated data, checkpoints, predictions, metrics, and figures are ignored by
Git. See [`docs/OUTPUTS.md`](docs/OUTPUTS.md) for the full artifact map.

## Using real data

Real inputs must follow the column contract in
[`docs/DATA_SCHEMA.md`](docs/DATA_SCHEMA.md). Most importantly:

1. harmonize alleles, genome build, chromosome labels, and effect directions;
2. create candidate-gene, tissue, and cell-type feature blocks in the documented
   order;
3. derive training labels independently of the features used to evaluate the
   model;
4. use chromosome- or locus-level holdouts to limit information leakage;
5. fit preprocessing on training data only. The saved checkpoint contains the
   training feature scaler and reuses it for validation, testing, inference,
   and ablation evaluation.

## Testing

```bash
python -m pip install -e ".[dev]"
python -m compileall -q sexreg_ad scripts run_pipeline.py
python -m pytest -q
```

## Citation and license

If you publish work based on this repository, cite the associated paper or
archive record once available. The code is released under the MIT License; see
[`LICENSE`](LICENSE).

## Contributing

Bug reports and focused pull requests are welcome. Please read
[`CONTRIBUTING.md`](CONTRIBUTING.md) before submitting changes.
