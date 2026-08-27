For more information on running FOCUS, please see the original Github repo's [wiki.](https://github.com/bogdanlab/focus/wiki)

Our paper uses a [slightly modified version of FOCUS](https://github.com/shugamoe/focus) with changes (some from myself, others from [Alvaro Barbeira](https://github.com/heroico)) to accomodate its use with
GTEx V8 models from predictdb.org.
# See https://github.com/shugamoe/focus for the specific FOCUS code changes.

## A Python 3.8 environment is required; this study uses FOCUS v0.7!

## Files 

|-- 000_run_format_gwas.pbs
|-- 000_run_format_gwas.py
   * These files format the relevant GWAS files from `../../00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/` into a format more suitable for FOCUS's built in "munging" (GWAS formatting) function.

|-- 001_run_patch_dbs_eqtl.pbs
|-- 001_run_patch_dbs_eqtl.py
|-- 002_run_patch_dbs_sqtl.pbs
|-- 002_run_patch_dbs_sqtl.py
  * This file was run to replace the rsID column with the varID column in the
  * predictdb.org databases, essentially a hack to accomodate the fact that we
  * had our LD data .bim files in ~ <chr_num>_<pos>_<a1>_<a2> sort of format.

|-- 010_run_focus_munge.pbs
  * Multiple versions of this file existed, depending on the case/control count
  * of the study. This is FOCUS's built in GWAS processing step.

|-- 020_run_focus_eqtl_sqtl_14tiss_one_by_one_db.sh
  * These files create the pre-requisite `.db` files needed to run FOCUS. A
  * single `.db` file is created for each (14) tissue for both expression and
  * splicing.

|-- 030_run_generate_by_tissue_runs.py
  * Generates runs of FOCUS, 1 run per chromosome per (14) tissues.

|-- 040_run_focus_eqtl_14tiss.sh
|-- 040_run_focus_sqtl_14tiss.sh
  * Submits the eqtl/sqtl based FOCUS runs in a way that doesn't blitz SLURM.

|-- 090_run_coalesce_by_tissue_output.py
  * Comboned FOCUS results.
  * Coallesce the individual outputs by chromosome/tissue into
  * a single file per study for convenience.



















