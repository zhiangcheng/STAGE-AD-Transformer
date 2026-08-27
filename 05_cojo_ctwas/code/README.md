# Overview

## COJO
# https://cloufield.github.io/GWASTutorial/18_Conditioning_analysis/
# gencode_v26_all.txt could be founded in hakyimlab.

010_run_cojo_eqtl.py
020_run_cojo_sqtl.py

  * Runs the COJO part of COJO/CTWAS, for all studies, expression and splicing.
  * Clumping R2 in use is `0.65`
  * We look 2MB left and right of a gene stat/end for previously reported GWAS variants.
  * We include ~10 or fewer variants we discovered in the METAL GWAS in the above GWAS variants.


## CTWAS
  * This folder contains conditional TWAS section of the COJO/CTWAS analysis. It relies on output from COJO. 
  
030_run_ctwas_prepare_eqtl.py
040_run_ctwas_prepare_sqtl.py
  
    * These files parse the COJO results, much like `00_prep_inputs_for_metaxcan/code/*parse*`
    * A key difference is that instead of having one result file per study, we have one result per gene of interest.


### CTWAS Splicing
  * If a folder is a "splicing" folder ("{study_name}")
  |-- 05_multi_tissue_intronxcan.jinja
  |-- 05_multi_tissue_intronxcan.sh
  |-- 05_multi_tissue_intronxcan.yaml
    * These files above are similar to those in `02_acat_eqtl_sqtl/code/029*`,
    * but they output results by gene/tissue list (either all or 11 tissues),
    * rather than by study/tissue list.

  |-- 051_combine_compare_mtiss_ixcan.py -> ../utils/051_combine_compare_mtiss_ixcan.py
    * This combines the above files into 1 file per study per tissue list.
    * Comparison is no longer done in this file, as it introduces some weird
    * circularity into the analysis that we can avoid.
  
### CTWAS Expression
  * If a folder is a "splicing" folder ("{study_name}_expression")
  |-- 06_1step_acat_and_combine.py -> ../utils/06_1step_acat_and_combine.py
    * This file calculates acat the same way as in
    * `02_acat_eqtl_sqtl/code/010_acat_eqtl*`, but is leaner code that doesn't
    * rely on MetaXCan code.