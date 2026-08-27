#!/bin/bash
### Preamble ###

directory=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/Males # This is where the plots will be output

p_threshold=5e-08 # p-value threshold for GWAS (where horizontal line will be placedon manhattan plot): 5e-08 for only common variants or 5e-09 for rare and common variants

files=${directory}/*_ad_male_sumstats_formatted_formetaanalysis.txt

#name=$(echo ${files} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/Males\///g' | sed 's/_ad_male_sumstats_formatted_formetaanalysis.txt//g')


### Run script ###

source activate GWAS
cd ${directory}


for FILE in ${files}
do
name=$(echo ${FILE} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/Males\///g' | sed 's/_ad_male_sumstats_formatted_formetaanalysis.txt//g')

# First change chr X for 23 so plots OK
awk '{ gsub("X", "23", $1) ; print }' ${FILE} > ${name}_temp.txt


Rscript --vanilla 04_Plots_Man_QQ.R ${name}_temp.txt ${name}_ad_male_sumstats_formatted_formetaanalysis ${p_threshold} 

rm ${name}_temp.txt

done
