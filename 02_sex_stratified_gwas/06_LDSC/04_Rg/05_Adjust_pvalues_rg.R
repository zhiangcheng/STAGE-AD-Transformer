# For rg results adjust p-value for multiple comparisons

directory = "/path/04_LDSC/SNPrg/"

rg_females = paste0(directory, "Females/rg_femaleMDD_correlation_table.txt", sep = "")

rg_males = paste0(directory, "Males/rg_maleMDD_correlation_table.txt", sep = "")


### LOAD PACKAGES ###
library(tidyverse)



### LOAD DATA ###

rg_females_df <- read.table(rg_females, header = T, stringsAsFactors = F)

rg_males_df <- read.table(rg_males, header = T, stringsAsFactors = F)


### Rg to previous MDD sumstats ###

rg_females_MDD <- rg_females_df %>% 
  filter (p2 == "AdamsMDD_sumstats" | p2 == "BloklandMDDfemale_sumstats" | p2 == "BloklandMDDmale_sumstats" | p2 == "SilveiraMDDfemale_sumstats" | p2 == "SilveiraMDDmale_sumstats") 

rg_males_MDD <- rg_males_df %>% 
  filter (p2 == "AdamsMDD_sumstats" | p2 == "BloklandMDDfemale_sumstats" | p2 == "BloklandMDDmale_sumstats" | p2 == "SilveiraMDDfemale_sumstats" | p2 == "SilveiraMDDmale_sumstats") 

rg_MDD <- bind_rows(rg_females_MDD, rg_males_MDD) %>% 
  mutate(padj = p.adjust(p, method = "BH", n = 10))


write.table(rg_MDD, paste0(directory, "rg_femaleMDD_and_maleMDD_with_MDD_correlation_table_padj.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

### Rg to other traits ###

rg_females_other <- rg_females_df %>% 
  filter (p2 != "AlsMDD_sumstats" & p2 != "AdamsMDD_sumstats" & p2 != "BloklandMDDfemale_sumstats" & p2 != "BloklandMDDmale_sumstats" & p2 != "SilveiraMDDfemale_sumstats" & p2 != "SilveiraMDDmale_sumstats" & p2 != "Allcohorts_male_MDD" & p2 != "BMI_female_sumstats" & p2 != "BMI_male_sumstats") 

rg_males_other <- rg_males_df %>% 
  filter (p2 != "AlsMDD_sumstats" & p2 != "AdamsMDD_sumstats" & p2 != "BloklandMDDfemale_sumstats" & p2 != "BloklandMDDmale_sumstats" & p2 != "SilveiraMDDfemale_sumstats" & p2 != "SilveiraMDDmale_sumstats" & p2 != "Allcohorts_female_MDD" & p2 != "BMI_female_sumstats" & p2 != "BMI_male_sumstats") 

rg_other <- bind_rows(rg_females_other, rg_males_other) %>% 
  mutate(padj = p.adjust(p, method = "BH", n = 22))


write.table(rg_other, paste0(directory, "rg_femaleMDD_and_maleMDD_with_other_traits_correlation_table_padj.txt"), col.names = T, row.names = F, quote = F, sep = "\t")


### Rg to other sex-specific BMI ###

rg_females_bmi <- rg_females_df %>% 
  filter (p2 == "BMI_female_sumstats" | p2 == "BMI_male_sumstats") 

rg_males_bmi <- rg_males_df %>% 
  filter (p2 == "BMI_female_sumstats" | p2 == "BMI_male_sumstats") 

rg_bmi <- bind_rows(rg_females_bmi, rg_males_bmi) %>% 
  mutate(padj = p.adjust(p, method = "BH", n = 4))


write.table(rg_bmi, paste0(directory, "rg_femaleMDD_and_maleMDD_with_sex_BMI_correlation_table_padj.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

