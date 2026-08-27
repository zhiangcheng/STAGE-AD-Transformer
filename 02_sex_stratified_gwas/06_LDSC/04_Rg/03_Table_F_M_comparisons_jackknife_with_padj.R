# Combine all rg comparisons into tables grouped by type of traits and adjust p-value

directory = "/path/04_LDSC/SNPrg/jack_knife/"


### LOAD PACKAGES ###
library(tidyverse)



### LOAD DATA ###
traits <- c(
  "AdamsMDD_sumstats",
  "BloklandMDDfemale_sumstats",
  "BloklandMDDmale_sumstats",
  "SilveiraMDDfemale_sumstats",
  "SilveiraMDDmale_sumstats",
  "ADHD_sumstats",
  "Anxiety_sumstats",
  "Bipolar_sumstats",
  "PTSD_sumstats",
  "Schizophrenia_sumstats",
  "BMI_sumstats",
  "Waist_to_hip_ratio_sumstats",
  "EA_sumstats",
  "Metabolic_syndrome_sumstats",
  "Drinks_per_week_sumstats",
  "Ever_smoked_regularly_sumstats",
  "BMI_female_sumstats",
  "BMI_male_sumstats"
)

rg_comparisons_df <- map(traits, function(trait) {
  read.table(paste0(directory, "Comparison_of_female_and_male_rg_with_", trait, ".txt"), header = TRUE, stringsAsFactors = FALSE)
}) %>%
  bind_rows()



### Rg to previous MDD sumstats ###

rg_comparisons_previous_MDD_df <- rg_comparisons_df %>%
  filter(trait %in% c(
    "AdamsMDD_sumstats",
    "BloklandMDDfemale_sumstats",
    "BloklandMDDmale_sumstats",
    "SilveiraMDDfemale_sumstats",
    "SilveiraMDDmale_sumstats"
  )) %>%
  mutate(padj = p.adjust(P, method = "BH", n = 5))

write.table(rg_comparisons_previous_MDD_df, paste0(directory, "Comparison_of_female_and_male_rg_with_previous_MDD.txt"), col.names = T, row.names = F, quote = F, sep = "\t")



### Rg to other traits ###

rg_comparisons_other_traits_df <- rg_comparisons_df %>%
  filter(trait %in% c(
    "ADHD_sumstats",
    "Anxiety_sumstats",
    "Bipolar_sumstats",
    "PTSD_sumstats",
    "Schizophrenia_sumstats",
    "BMI_sumstats",
    "Waist_to_hip_ratio_sumstats",
    "EA_sumstats",
    "Metabolic_syndrome_sumstats",
    "Drinks_per_week_sumstats",
    "Ever_smoked_regularly_sumstats"
  )) %>%
  mutate(padj = p.adjust(P, method = "BH", n = 11))

write.table(rg_comparisons_other_traits_df, paste0(directory, "Comparison_of_female_and_male_rg_with_other_traits.txt"), col.names = T, row.names = F, quote = F, sep = "\t")



### Rg to sex-specific BMI ###

rg_comparisons_sex_specific_bmi_df <- rg_comparisons_df %>%
  filter(trait %in% c(
    "BMI_female_sumstats",
    "BMI_male_sumstats"
  )) %>%
  mutate(padj = p.adjust(P, method = "BH", n = 2))

write.table(rg_comparisons_sex_specific_bmi_df, paste0(directory, "Comparison_of_female_and_male_rg_with_sex_specific_bmi.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

