# Identify putative causal risk loci for identified regions by gwas-pw


directory="/path/06_GWAS_pw/femaleMDD_maleMDD/"

femaleMDD_maleMDD_SNP <- "/path/06_GWAS_pw/femaleMDD_maleMDD/GWAS_pw_femaleMDD_maleMDD.bfs.gz"

femaleMDD_maleMDD_region <- "/path/06_GWAS_pw/femaleMDD_maleMDD/GWAS_pw_femaleMDD_maleMDD.segbfs.gz"





### LOAD PACKAGES ###
library(tidyverse)

### READ IN DATA ###

femaleMDD_maleMDD_SNP_df <- read.table(gzfile(femaleMDD_maleMDD_SNP), header = T, stringsAsFactors = F)

femaleMDD_maleMDD_region_df <- read.table(gzfile(femaleMDD_maleMDD_region), header = T, stringsAsFactors = F)




### CAUSAL VARIANTS FOR SHARED REGIONS ###

# Within each chunk identified as having a shared causal variant between the 2 traits, identify which is the likely causal variant (SNP with the highest PPA_3 = posterior probability that this SNP is the causal one under model 3)

shared_regions <- femaleMDD_maleMDD_region_df %>% 
  filter(PPA_3 > 0.5)

# Top SNP per shared region
Top_SNP_per_chunk_shared <- femaleMDD_maleMDD_SNP_df %>%
  filter(chunk %in% shared_regions$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_3) %>%
  slice_head(n = 1) %>%
  ungroup()

outfile <- paste(directory, "Top_SNP_per_shared_region.txt", sep="")
write.table(Top_SNP_per_chunk_shared, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# SNPs per shared region with PPA_3 > 0.5
PPA_05_SNP_per_chunk_shared <- femaleMDD_maleMDD_SNP_df %>%
  filter(chunk %in% shared_regions$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_3 > 0.5) %>%
  ungroup()

outfile <- paste(directory, "SNP_with_PPA_above_05_per_shared_region.txt", sep="")
write.table(PPA_05_SNP_per_chunk_shared, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



### CAUSAL VARIANTS FOR FEMALE ONLY REGIONS ###

# Within each chunk identified as having a causal variant in females and not males, identify which is the likely causal variant (SNP with the highest PPA_1 = posterior probability that this SNP is the causal one under model 3)

female_regions <- femaleMDD_maleMDD_region_df %>% 
  filter(PPA_1 > 0.5)

# Top SNP per shared region
Top_SNP_per_chunk_female <- femaleMDD_maleMDD_SNP_df %>%
  filter(chunk %in% female_regions$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_1) %>%
  slice_head(n = 1) %>%
  ungroup()

outfile <- paste(directory, "Top_SNP_per_female_region.txt", sep="")
write.table(Top_SNP_per_chunk_female, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# SNPs per female region with PPA_1 > 0.5
PPA_05_SNP_per_chunk_female <- femaleMDD_maleMDD_SNP_df %>%
  filter(chunk %in% female_regions$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_1 > 0.5) %>%
  ungroup()

outfile <- paste(directory, "SNP_with_PPA_above_05_per_female_region.txt", sep="")
write.table(PPA_05_SNP_per_chunk_female, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


### CAUSAL VARIANTS FOR MALE ONLY REGIONS ###
male_regions <- femaleMDD_maleMDD_region_df %>% 
  filter(PPA_2 > 0.5)
# No male only regions


### CAUSAL VARIANTS FOR SHARED REGIONS WITH DIFFERENT CAUSAL VARIANT ###
shared_diff_regions <- femaleMDD_maleMDD_region_df %>% 
  filter(PPA_4 > 0.5)
# None
