#!/bin/bash
### Preamble ###

directory=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/

files=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/GxS/nodc/*_ad_GxS_nodc_sumstats_formatted.txt

#name=$(echo ${files} | sed ''s/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\///g' | sed 's/_ad_GxS_nodc_sumstats_formatted.txt//g')


### Run script ###

source activate GWAS
cd ${directory}

for FILE in ${files}
do

name=$(echo ${FILE} | sed 's/\/media\/desk15\/iy2120\/Project2026\/Myproject_AD\/02_sex_stratified_gwas\/010_Format_sumstats\/GxS\/nodc\///g' | sed 's/_ad_GxS_nodc_sumstats_formatted.txt//g')

echo "Checking file: ${FILE}" > ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt


## Check how many SNPs did not match to marker name ##

echo "Number rows missing (including non-SNPs) - no match in dbSNP ref file:" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
grep -c "missing" ${FILE} >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt

echo "Number rows missing (SNPs only) - no match in dbSNP ref file:" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
awk '$16 == "missing" && length($3) == 1 && length($4) == 1' ${FILE} | wc -l >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt

echo "Number rows with \"no_match\" (SNP present in the reference dbSNP file but not in the particular build using):" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
grep -c "no_match" ${FILE} >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt

echo "Number rows with \"no_rsid\" (SNP present in one build, mapped to the other build and coords found but no matching rsid found):" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
grep -c "no_rsid" ${FILE} >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt

echo "No. SNPs in final file after matching, formatting and transforming:" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
wc -l ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt




## Check min and max of each column to see if anything weird ##

# Before removing missing lines

echo "Min and max of each column before removing missing lines" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
awk 'NR==1 {
        for (i=1; i<=NF; i++)
            headers[i] = $i  # Store column headers in an array
     }
     NR > 1 {
        for (i=1; i<=NF; i++) {
            if (NR == 2 || $i < min[i]) {  # Update min value for each column
                min[i] = $i
                minHeader[i] = headers[i]  # Store corresponding header
            }
            if (NR == 2 || $i > max[i]) {  # Update max value for each column
                max[i] = $i
                maxHeader[i] = headers[i]  # Store corresponding header
            }
        }
     }
     END {
        for (i=1; i<=NF; i++) {
            print "Min value for column " headers[i] ": " min[i]
            print "Max value for column " headers[i] ": " max[i]
        }
     }' ${FILE} >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt


# After removing missing lines
echo "Min and max of each column after removing missing lines" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
awk 'NR==1 {
        for (i=1; i<=NF; i++)
            headers[i] = $i  # Store column headers in an array
     }
     NR > 1 {
        for (i=1; i<=NF; i++) {
            if (NR == 2 || $i < min[i]) {  # Update min value for each column
                min[i] = $i
                minHeader[i] = headers[i]  # Store corresponding header
            }
            if (NR == 2 || $i > max[i]) {  # Update max value for each column
                max[i] = $i
                maxHeader[i] = headers[i]  # Store corresponding header
            }
        }
     }
     END {
        for (i=1; i<=NF; i++) {
            print "Min value for column " headers[i] ": " min[i]
            print "Max value for column " headers[i] ": " max[i]
        }
     }' ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis.txt >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt


# After transforming beta and SE to logistic scale
echo "Min and max of each column after removing missing lines and transforming beta and SE to logistic scale" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
awk 'NR==1 {
        for (i=1; i<=NF; i++)
            headers[i] = $i  # Store column headers in an array
     }
     NR > 1 {
        for (i=1; i<=NF; i++) {
            if (NR == 2 || $i < min[i]) {  # Update min value for each column
                min[i] = $i
                minHeader[i] = headers[i]  # Store corresponding header
            }
            if (NR == 2 || $i > max[i]) {  # Update max value for each column
                max[i] = $i
                maxHeader[i] = headers[i]  # Store corresponding header
            }
        }
     }
     END {
        for (i=1; i<=NF; i++) {
            print "Min value for column " headers[i] ": " min[i]
            print "Max value for column " headers[i] ": " max[i]
        }
     }' ${name}_ad_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt
echo "" >> ${name}_ad_GxS_nodc_sumstats_check_formatted_file.txt



# Check lowest p-value retained after matching consistent marker
# Check QC on MAF and info done correctly

done
