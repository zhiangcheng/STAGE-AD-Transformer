# Input data schema

The training and inference entry points read a Parquet table with one row per
variant. Column order is not significant, but names and token indices are.
`scripts/01_simulate_data.py` is the executable schema example.

## Metadata and labels

| Column | Type | Meaning | Needed for inference? |
|---|---|---|---|
| `variant_id` | string | stable variant identifier | yes |
| `rsid` | string | dbSNP-style identifier | yes |
| `CHR` | string | `1`–`22` or `X` | yes |
| `BP` | integer | base-pair position in one declared genome build | yes |
| `nearest_gene` | string | nearest annotated gene | yes |
| `sex_label` | string | human-readable five-class label | current predictor output metadata |
| `sex_label_id` | integer | class ID 0–4 | evaluation/training |
| `is_positive` | 0/1 | AD-risk binary target | evaluation/training |
| `sex_context_id` | integer | 0 male, 1 female, 2 combined/neutral | yes |
| `target_gene_label` | integer | correct candidate index 0–5 | evaluation/training |
| `top_tissue_label` | integer | correct tissue index 0–8 | evaluation/training |
| `top_cell_label` | integer | correct cell index 0–6 | evaluation/training |
| `mechanism_label` | integer | correct mechanism index 0–7 | evaluation/training |
| `chr_type` | string | autosome, PAR, or X_nonPAR | prediction metadata |
| `xci_status` | string | autosome, PAR, escape, subject, or variable | prediction metadata |

The current inference script expects label columns because the same dataset
class supports evaluation. For unlabeled deployment data, add sentinel labels
within their valid ranges or extend the dataset with a dedicated inference
mode before use. Do not interpret metrics computed from sentinel labels.

## Variant feature block

The exact ordered list is `sexreg_ad.constants.VARIANT_FEATURES`:

- quality and LD: `MAF`, `INFO`, `LD_score`, `distance_to_tss_kb`;
- association: `Z_*`, `neglog10P_*`, `Beta_*`, sex-effect differences;
- fine-mapping: `PIP_total`, `PIP_female`, `PIP_male`, `PIP_interaction`;
- colocalization: `coloc_brain_pph4`, `coloc_blood_pph4`,
  `coloc_immune_pph4`;
- functional evidence: `deepsea_delta`, `enformer_brain_delta`,
  `sei_regulatory_score`, `cadd_score`, `conservation_score`;
- chromosome X/XCI indicators: `x_is_chrX`, `x_is_nonPAR`, `xci_escape`,
  `xci_subject`, `xci_variable`.

All variant features are required. Missing or misspelled columns raise a clear
schema error before model execution.

## Candidate-gene token block

There are six candidate genes, indexed `0`–`5`. Each candidate needs
`gene{k}_name` plus these numeric fields:

`distance_kb`, `eqtl_z`, `sqtl_z`, `pqtl_z`, `twas_z`, `coloc_pph4`,
`brain_expression`, `ad_de_z`, `constraint_score`.

Example: `gene0_name`, `gene0_distance_kb`, ..., `gene5_constraint_score`.

## Tissue token block

There are nine tissues in the fixed order defined by `TISSUES`. For each index
`t`, provide:

`tissue{t}_eqtl_z`, `tissue{t}_sqtl_z`, `tissue{t}_coloc_pph4`,
`tissue{t}_enhancer`, `tissue{t}_ldsc_enrichment`,
`tissue{t}_expression`.

## Cell-type token block

There are seven cell types in the fixed order defined by `CELL_TYPES`. For each
index `c`, provide:

`cell{c}_marker_score`, `cell{c}_ad_de_score`,
`cell{c}_enhancer_score`, `cell{c}_expression_specificity`.

## Preprocessing contract

- Convert infinities and missing numeric values to zero before scaling.
- Fit each feature's mean and standard deviation on the training split only.
- Save the scaler with the model checkpoint.
- Reuse exactly the saved scaler for validation, test, prediction, and ablation.
- Keep feature definitions and token order identical between training and use.

The repository implements this contract in `SexRegV3Dataset` and checkpoint
creation. Checkpoints from the older code without a `scaler` field remain
loadable, but inference then falls back to fitting the scaler on the input and
prints a warning.
