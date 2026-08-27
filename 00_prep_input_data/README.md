# Preparing AD GWAS Summary Statistics for Subsequent Analysis

## Scope and manuscript alignment

The manuscript's primary discovery analysis is a sex-aware genome-wide
association study of Alzheimer disease (AD), not the three-study,
sex-combined meta-analysis encoded by the legacy scripts in this directory.
The primary discovery dataset comprises 449,335 participants of genetically
inferred European ancestry from five cohorts:

- Alzheimer's Disease Genetics Consortium (ADGC);
- Alzheimer's Disease Sequencing Project (ADSP);
- European Alzheimer & Dementia Biobank (EADB);
- UK Biobank (UKBB); and
- All of Us.

The manuscript reports the following primary discovery counts:

| Analysis stratum | Cases | Controls | Total participants |
| --- | ---: | ---: | ---: |
| Female | 109,525 | 149,773 | 259,298 |
| Male | 69,423 | 120,614 | 190,037 |
| Combined discovery sample | 178,948 | 270,387 | 449,335 |

The corresponding sex-specific effective sample sizes are 109,349 for females
and 82,661 for males. The manuscript also describes independent replication in
230,987 participants from non-overlapping European- and East Asian-ancestry
cohorts.

Four frozen GWAS summary-statistic contexts are required for the manuscript's
cross-tissue analyses:

1. overall/sex-combined AD;
2. female-specific AD;
3. male-specific AD; and
4. SNP-by-sex interaction (GS interaction).

This directory does **not** generate the cohort-level fastGWA results or the
primary five-cohort sex-stratified and GS-interaction meta-analyses. Those
summary statistics must be generated and quality-controlled upstream. This
directory contains a separate legacy workflow that combines UKB, Kunkle 2019,
and PGC-ALZ 2021 data to produce an autosomal, sex-combined MetaXcan input.
That workflow can be retained as an auxiliary overall-reference analysis, but
it must not be described as the manuscript's primary discovery GWAS.

## Manuscript analysis flow

The manuscript-aligned downstream workflow is:

```text
Five-cohort AD GWAS
|-- overall/sex-combined summary statistics
|-- female-specific summary statistics
|-- male-specific summary statistics
`-- SNP-by-sex interaction summary statistics
        |
        v
Allele and genome-build harmonization (GRCh38)
        |
        v
Summary imputation to variants represented in the prediction models
        |
        |-- expression-based S-PrediXcan using eQTL models
        `-- splicing-based S-PrediXcan using sQTL models
        |
        v
Cross-tissue ACAT combination
        |
        |-- conditional TWAS after COJO adjustment
        `-- FOCUS gene-level fine-mapping
```

All genomic coordinates are harmonized to GRCh38. Effect alleles must be
aligned consistently across GWAS, QTL, LD-reference, and prediction-model
files before any downstream analysis.

## Upstream GWAS requirements from the manuscript

Sex-stratified GWAS are conducted separately in females and males using
fastGWA mixed models. GS-interaction models include genotype, sex, and a
genotype-by-sex interaction term. Autosomes and chromosome X are analyzed,
with the non-pseudoautosomal region of chromosome X evaluated under both no
and full dosage-compensation assumptions.

The primary GWAS filters retain variants with:

- minor allele frequency greater than 0.01; and
- imputation quality greater than 0.6.

Cohort-level summary statistics are combined by inverse-variance fixed-effect
meta-analysis. The manuscript's implementation uses standard-error-weighted
METAL without genomic control and retains variants present in at least three
of the five discovery cohorts. Random-effects models, Cochran's Q, I-squared,
and leave-one-cohort-out analyses are used to assess heterogeneous loci.

For the GS analysis, linear mixed-model estimates are converted to the
log-odds-ratio scale before meta-analysis, and the corresponding standard error
is obtained from the converted effect and its Z score. These operations belong
to the upstream GWAS workflow and are not implemented by the scripts in this
directory.

## Cross-tissue analyses described in the manuscript

The manuscript uses two joint-tissue TWAS strategies across 19 tissue types:

- **Expression-based TWAS:** S-PrediXcan is run for each gene in each tissue
  using eQTL prediction models, followed by ACAT across tissues.
- **Splicing-based TWAS:** S-PrediXcan is run for sQTL models; signals are
  combined first within intron-excision clusters and then across tissues using
  ACAT.

The primary ACAT analysis weights the brain tissues collectively as much as
all non-brain tissues combined. Equal weighting across all tissues is used as
a sensitivity analysis. Both strategies are applied to the overall,
female-specific, male-specific, and GS-interaction GWAS summary statistics.

For conditional TWAS, eQTL and sQTL effect estimates are conditioned with COJO
on genome-wide significant AD index variants located within 2 Mb of a gene's
transcription start or stop site. FOCUS is then run separately for expression
and splicing models in each tissue. The manuscript defines 90% credible gene
sets and prioritizes candidate causal genes with a posterior inclusion
probability of at least 0.9.

> **Tissue-count warning:** some downstream scripts and file names in this code
> snapshot still refer to `14tiss`, whereas the current manuscript specifies
> 19 tissue types. Freeze and document the final 19-tissue list before public
> release, and update the downstream lists and loops accordingly.

## Legacy auxiliary workflow in this directory

The files in this directory implement the following sex-combined, autosomal
workflow:

1. meta-analyze UKB, Kunkle 2019, and PGC-ALZ 2021 with METAL;
2. harmonize the METAL output and lift coordinates from hg19 to GRCh38;
3. impute autosomal variants against a European 1000 Genomes reference panel;
4. merge the observed and imputed regional results into one MetaXcan-ready
   file; and
5. optionally add rsIDs.

The harmonization and imputation commands follow the
[Summary-GWAS-Imputation workflow](https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation).

### Legacy input files

| Study | Input file | Effect | Standard error | Effect/non-effect alleles | Frequency |
| --- | --- | --- | --- | --- | --- |
| UK Biobank AD | `UKB_AD_clean.txt.gz` | `REGENIE_BETA` | `REGENIE_SE` | `Alt` / `Ref` | `MAF` |
| Kunkle 2019 Stage 1 | `Kunkle_Stage1_gwas_clean.txt.gz` | `Beta` | `SE` | `Effect_allele` / `Non_Effect_allele` | `eaf` |
| PGC-ALZ 2021 | `PGCALZ_gwas_clean.txt.gz` | `Beta` | `SE` | `Effect_allele` / `Non_Effect_allele` | `eaf` |

All three files must contain a consistently constructed `MarkerName` field.
The METAL configuration uses a standard-error-weighted model, reports allele
frequency ranges, performs heterogeneity analysis, and does not apply genomic
control.

The values currently hard-coded in
`020_parse_final_metal_ad_kunkle_pgcalz_ukb.pbs` are 47,525 cases and
1,121,785 controls. These values belong to the legacy auxiliary meta-analysis,
are not the manuscript's discovery sample sizes, and must be verified against
the frozen source-data ledger and overlap audit before reuse.

### Step 1: legacy METAL meta-analysis

Run:

```bash
bash 010_run_metal_ad_kunkle_pgcalz_ukb.pbs
```

This script uses `metal_ad_kunkle_pgcalz_ukb.txt` and writes:

```text
output/metal/final_metal_ad_kunkle_pgcalz_ukb1.txt
```

### Step 2: harmonization and liftover

Run:

```bash
bash 020_parse_final_metal_ad_kunkle_pgcalz_ukb.pbs
```

The script maps METAL columns to the Summary-GWAS-Imputation schema, splits
`MarkerName` into chromosome, position, and allele fields, lifts coordinates
from hg19 to GRCh38, and aligns variants to GTEx v8 European reference
metadata. It writes:

```text
output/metaxcan_gwas_parse/final_metal_ad_kunkle_pgcalz_ukb1_EUR.txt.gz
```

`test_parse_pgcalz_gwas_clean.pbs` is an optional diagnostic for PGC-ALZ
parsing and is not part of the main auxiliary workflow.

### Step 3: autosomal summary imputation

Copy the marked Bash block in `030_run_summary_imputation.pbs` into
`030_run_summary_imputation.sh`, then run the PBS wrapper. The current script
processes chromosomes 1-22 in 10 sub-batches per chromosome with:

- a 100-kb window;
- `regularization = 0.1`;
- `frequency_filter = 0.01`;
- `parsimony = 9`; and
- standardized dosages.

Regional results are written to:

```text
output/metaxcan_gwas_impute_regions/
```

This loop is autosomal only. It does not reproduce the manuscript's separate
chromosome-X association analyses.

### Step 4: post-processing

Run:

```bash
bash 040_post_process_final_metal_ad_kunkle_pgcalz_ukb.pbs
```

This merges the harmonized observed statistics with all regional imputation
results and writes:

```text
output/metaxcan_gwas_imputed/final_metal_ad_kunkle_pgcalz_ukb.txt.gz
```

### Step 5: optional rsID annotation

`050_add_rsid_col_to_final_metal_ad_kunkle_pgcalz_ukb_results.pbs` launches the
corresponding R script. The R file still contains breast-cancer paths and
lookup schemas and is therefore retained only as a legacy template. Replace
`YOUR_METAL_FILE`, `OUT_METAL_PATH`, and all rsID reference inputs before use.

## Adapting this workflow to the manuscript inputs

For a manuscript-aligned run, start from the four frozen upstream summary
statistics rather than rerunning the legacy three-study METAL step. For each
analysis context:

1. verify the phenotype, cohort composition, genome build, sample size,
   effect-allele convention, and overlap exclusions against the provenance
   ledger;
2. map columns to the harmonized schema used by
   Summary-GWAS-Imputation;
3. set the correct stratum-specific sample-size fields rather than copying the
   legacy hard-coded values;
4. lift to GRCh38 only when the frozen input is not already in GRCh38;
5. align alleles against the same European reference metadata and prediction
   models used downstream;
6. run summary imputation separately for overall, female, male, and
   GS-interaction statistics; and
7. freeze the output name, checksum, input version, and software version before
   S-PrediXcan, ACAT, COJO, or FOCUS analyses.

Do not combine the female, male, and interaction files in METAL at this stage;
they represent distinct analysis contexts and must remain separate throughout
TWAS and validation.

## Software and reference data

The auxiliary workflow requires:

- METAL;
- Python 3 in the `TWAS` Conda environment;
- `gwas_parsing.py`, `gwas_summary_imputation.py`, and
  `gwas_summary_imputation_postprocess.py` from Summary-GWAS-Imputation;
- `hg19ToHg38.over.chain.gz`;
- `gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz`;
- `eur_ld.bed.gz`;
- European 1000 Genomes genotype parquet files and
  `variant_metadata.parquet`; and
- R with `tidyverse`, `data.table`, and `glue` for the optional rsID step.

The downstream manuscript analyses additionally require S-PrediXcan, the
frozen eQTL and sQTL prediction models, ACAT, GCTA-COJO, tissue-specific LD
matrices, and FOCUS.

## Path configuration

The scripts use absolute paths under:

```text
$HOME/Project2026/Myproject_AD/00_prep_inputs_for_metaxcan/
```

This repository snapshot stores the files in `00_prep_input_data/`. Reproduce
the expected runtime layout or replace the hard-coded project paths before
running the workflow. Some shared software and reference paths still point to
`$HOME/TWAS/Breast-cancer-Example/`; these paths must also be replaced or
documented for a public release.

Expected runtime subdirectories are:

```text
00_prep_inputs_for_metaxcan/
|-- code/
|-- input/
|-- output/
|   |-- metal/
|   |-- metaxcan_gwas_parse/
|   |-- metaxcan_gwas_impute_regions/
|   `-- metaxcan_gwas_imputed/
`-- temp/
```

## Reproducibility checks before release

- Confirm that the README sample counts match the frozen manuscript and
  Supplementary Data.
- Resolve any conflicting sample counts elsewhere in the manuscript before
  assigning a release tag.
- Document all cohort-overlap exclusions.
- Record whether each input is overall, female, male, or GS interaction.
- Record genome build, effect allele, non-effect allele, frequency definition,
  sample size, case/control counts, and imputation-quality field.
- Confirm that the manuscript's 19-tissue list matches the tissue lists used by
  S-PrediXcan and ACAT.
- Keep chromosome-X results separate from this autosomal summary-imputation
  workflow.
- Verify that no legacy breast-cancer path or phenotype label remains in an
  executable manuscript workflow.
