#!/bin/bash
### Preamble ###

directory=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/ # This is where the plots will be output

dbSNP=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/sex_stratified_data/dbSNP155/NCBI_Hsapiens_dbSNP155_GRCh37_GRCh38.p13_split_multiallelics_MAF0_0001_ALL_matches.txt

files=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/*_ad_GxS_nodc_sumstats.txt

#name=$(echo ${files} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\///g' | sed 's/_build37_ad_GxS_nodc_sumstats.txt//g')

### Run script ###

source activate GWAS
cd ${directory}

# Check if the input file name contains "build37" or "build38" and run the appropriate block of code
# MarkerName is then created based on the dbSNP ref file so it is consistent across all cohort sumstats so corrcet matching occurs in the meta-analysis
# Note that MarkerName is in the format CHR:BP:Allele:Allele but doesn't necessarily match to A1:A2 (alleles are not flipped - but this is not a problem as METAL deals with allele flipping)
for FILE in ${files}
do
name=$(echo ${FILE} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\///g' | sed 's/_build37_ad_GxS_nodc_sumstats.txt//g')

if [[ "${FILE}" == *"build37"* ]]; then
    echo "Processing file for build 37: ${FILE}"
    awk 'NR==FNR {
           marker[$1,$2,$4,$5]=$10;
            rsid[$1,$2,$4,$5]=$3;
            next
         }
         {
            if (FNR == 1)
                print $0,"MARKER_build37","rsID_build37";
            else if (($1,$2,$3,$4) in marker)
                print $0,marker[$1,$2,$3,$4],rsid[$1,$2,$3,$4];
            else if (($1,$2,$4,$3) in marker)
                print $0,marker[$1,$2,$4,$3],rsid[$1,$2,$4,$3];
            else
                print $0,"missing","missing";
         }' ${dbSNP} ${FILE} > ${name}_ad_GxS_nodc_sumstats_formatted.txt


elif [[ "${FILE}" == *"build38"* ]]; then
    echo "Processing file for build 38: ${FILE}"
    awk 'NR==FNR {
            marker[$6,$7,$8,$9]=$10;
            rsid[$6,$7,$8,$9]=$3;
            next
         }
         {
            if (FNR == 1)
                print $0,"MARKER_build37","rsID_build37";
            else if (($1,$2,$3,$4) in marker)
                print $0,marker[$1,$2,$3,$4],rsid[$1,$2,$3,$4];
            else if (($1,$2,$4,$3) in marker)
                print $0,marker[$1,$2,$4,$3],rsid[$1,$2,$4,$3];
            else
                print $0,"missing","missing";
         }' ${dbSNP} ${FILE} > ${name}_ad_GxS_nodc_sumstats_formatted.txt
else
    echo "Error: Unrecognized build in file name: ${FILE}. Skipping."
fi


## Remove lines with no match to dbSNP reference file ##

awk '($16 != "missing" && !($16 ~ /^no_match/) && $17 != "no_rsid") {print}' ${name}_ad_GxS_nodc_sumstats_formatted.txt > ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis.txt

done
