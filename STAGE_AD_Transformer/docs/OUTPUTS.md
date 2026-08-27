# Output reference

Generated files are deliberately excluded from Git. Archive large public
artifacts in a suitable data repository and link them from a release or paper.

## Data and splits

| Path | Description |
|---|---|
| `data/processed/simulated_variants.parquet` | complete synthetic table |
| `data/processed/simulated_by_chr/` | chromosome-partitioned synthetic data |
| `data/splits/train.parquet` | training split |
| `data/splits/valid.parquet` | validation split |
| `data/splits/test.parquet` | held-out chromosome test split |
| `data/splits/split_summary.csv` | sizes, class balance, and chromosomes |

## Checkpoints and metrics

| Path | Description |
|---|---|
| `results/checkpoints/sexreg_ad_pretrained.pt` | pretrained weights, arguments, scaler |
| `results/checkpoints/sexreg_ad_finetuned.pt` | selected fine-tuned model, arguments, scaler |
| `results/metrics/pretraining_history.csv` | epoch-level pretraining metrics |
| `results/metrics/finetuning_history.csv` | epoch-level fine-tuning metrics |
| `results/metrics/transformer_test_metrics.csv` | final held-out metrics |
| `results/metrics/baseline_metrics.csv` | traditional-method comparison |
| `results/metrics/ablation_metrics.csv` | inference-time ablation results |

## Predictions

`results/predictions/predictions.csv` and `.parquet` contain variant metadata,
AD-risk probability, five-class probabilities, predicted target-gene index and
name, predicted tissue, cell type, mechanism, and female/male context scores.

## Final tables

The directory `results/final_outputs/` contains:

1. `00_model_registry.csv`
2. `01_variant_ad_risk_priority_table.csv`
3. `02_sex_class_prediction_table.csv`
4. `03_top_target_genes_per_locus.csv`
5. `04_top_tissues_per_locus.csv`
6. `05_top_cell_types_per_locus.csv`
7. `06_model_attribution_by_domain.csv`
8. `07_benchmark_against_traditional_methods.csv`
9. `08_top_ranked_novel_loci.csv`
10. `09_aim3_candidate_mechanisms.csv`
11. `MANIFEST_final_deliverables.csv`

## Figures

`results/figures/` contains matching PNG and PDF versions of ten figure groups:
four Manhattan plots, model benchmark, ablation analysis, sex-context scatter,
cell-type heatmap, chromosome-X lollipop plot, and regulatory network.
