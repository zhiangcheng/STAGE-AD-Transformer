# Analysis plan

This document records the intended scientific analysis. It is not evidence
that the synthetic demonstration validates biological or clinical claims.

## Prediction tasks

1. Binary AD-risk classification.
2. Five-class sex-effect classification: null, shared, female-biased,
   male-biased, and sex-interaction.
3. Candidate target-gene ranking.
4. Tissue ranking.
5. Cell-type ranking.
6. Regulatory-mechanism classification.

## Architecture

- CLS and sex-context tokens.
- One learned token per scalar variant feature.
- Six candidate-gene tokens.
- Nine tissue tokens and seven cell-type tokens.
- Token-type embeddings and a Transformer encoder.
- Gene, tissue, and cell cross-attention into the CLS representation.
- Six task heads plus masked variant-feature reconstruction.

## Training and evaluation protocol

- Fit feature scaling on the training split only and save it in checkpoints.
- Hold out complete chromosomes, including chromosome X, for testing.
- Select the fine-tuned model using validation AUPRC.
- Report AUROC, AUPRC, F1, balanced accuracy, Brier score, ECE, and calibration
  slope for the binary task.
- Report macro/weighted F1 and per-class AUPRC for the five-class task.
- Report top-1/top-3 accuracy and MRR for ranking tasks.
- Compare fixed-score and learned baselines on the same held-out test data.
- Label inference-time feature masking as ablation, not retraining ablation.

## Real-data requirements

Before biological interpretation, replace all synthetic labels and features
with harmonized GWAS, QTL/TWAS, fine-mapping, functional-sequence,
single-cell, XCI, and external validation evidence. Record cohort inclusion,
genome build, allele harmonization, ancestry, covariates, label construction,
and all leakage controls in the associated study protocol.
