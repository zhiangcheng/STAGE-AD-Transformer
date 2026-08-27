#!/bin/bash
### Preamble ###

# Edit below variables for your specific LDSC

directory=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/ # working directory

reference=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/sex_stratified_data/LDSC_reference/ # location of reference directory

files=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/Munged_sumstats/*_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged.sumstats.gz

#name=$(echo ${files} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\/Munged_sumstats\///g' | sed 's/_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged.sumstats.gz//g')


### Submit script ###

source activate LDSC
cd ${directory}

for FILE in ${files}
do

name=$(echo ${FILE} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\/Munged_sumstats\///g' | sed 's/_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged.sumstats.gz//g')
echo ${name}

ldsc.py \
 --h2 ${FILE} \
 --ref-ld-chr ${reference}eur_w_ld_chr/ \
 --w-ld-chr ${reference}eur_w_ld_chr/ \
 --out ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_h2

# eur_w_ld_chr directory provided by the software

# Create output with info on intercept

echo "LDSC Intercept and SD:" >> ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSCintercept.txt
grep "Intercept:" ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_h2.log >> ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSCintercept.txt

intercept=$(grep "Intercept:" ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_h2.log | awk '{print $2}')
sd=$(grep "Intercept:" ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_h2.log | awk '{print $3}' | sed 's/[()]//g')
lower_ci=$(echo "$intercept - (1.96 * $sd)" | bc)
upper_ci=$(echo "$intercept + (1.96 * $sd)" | bc)
echo "LDSC Intercept 95% CI: Lower = $lower_ci, Upper = $upper_ci" >> ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSCintercept.txt

# In *_LDSC_h2.log file the intercept and the SD (in brackets) is reported
# Calculate the intercept's 95% CI
#  95% CI = 1.96 * SD
#  Therefore, intercept +/- (SD x 1.96) gives 95% CI
# Intercept not sig dif from 1 suggests lambda inflation is true polygenicity and not due to confounding.
# Intercept sig > 1 suggests confounding
# Intercept sig < 1 suggests overcorrection

done
