### Packages ###

library(tidyverse)

### Read In data ###

females <- read.table(file = "/path/04_LDSC/SNPrg/Females/rg_femaleMDD_correlation_table.txt", header = T)

males <- read.table(file = "/path/04_LDSC/SNPrg/Males/rg_maleMDD_correlation_table.txt", header = T)

rg_df <- females %>% 
  full_join(males) %>% 
  filter(p1 != "Allcohorts_male_MDD" | p2 != "Allcohorts_female_MDD")

### Test if the rg between males and females is sig different than 1 ###

rg_female_male_df <- rg_df %>% 
  filter(p1 == "Allcohorts_female_MDD" & p2 == "Allcohorts_male_MDD")

rg_female_male <- rg_female_male_df$rg

SE_female_male <- rg_female_male_df$se

# Z-score
Z <- (1 - rg_female_male) / SE_female_male

# P-value from a standard normal distribution that rg female - male is significantly below 1
p <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = F)

# Print results to file
cat("Statistical test that the rg between female and male MDD GWAS is significantly below 1: Z-score =", Z, "p =", p, "\n", file = "Rg_female_male_different_to_1.txt")


### Test if rg between female MDD and a trait is significantly different to male MDD and the same trait ###

my_function <- function(trait) {
  
  female_rg_df <- rg_df %>% 
    filter(p1 == "Allcohorts_female_MDD" & p2 == trait)
  
  female_rg <- female_rg_df$rg
  female_se <- female_rg_df$se
  
  male_rg_df <- rg_df %>% 
    filter(p1 == "Allcohorts_male_MDD" & p2 == trait)
  
  male_rg <- male_rg_df$rg
  male_se <- male_rg_df$se
  
  Z <- (female_rg - male_rg) / sqrt((female_se)^2 + (male_se)^2)
  
  p <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = F)
  
  result <- tibble(
    trait = trait,
    Z = Z,
    p = p
  )
}


trait_list_MDD <- c("AdamsMDD_sumstats", "BloklandMDDfemale_sumstats", "BloklandMDDmale_sumstats", "SilveiraMDDfemale_sumstats", "SilveiraMDDmale_sumstats")

results_df_MDD <- map_dfr(trait_list_MDD, my_function)

results_df_MDD <- results_df_MDD %>% 
  mutate(P_adjusted = p.adjust(p, method = "BH", n = 5)) 

write_delim(results_df_MDD, file = "Compare_female_rg_to_male_rg_for_MDD.txt")




trait_list_other <- c("ADHD_sumstats", "Anxiety_sumstats", "Bipolar_sumstats", "PTSD_sumstats", "Schizophrenia_sumstats", "BMI_sumstats", "Waist_to_hip_ratio_sumstats", "EA_sumstats", "Metabolic_syndrome_sumstats","Drinks_per_week_sumstats","Ever_smoked_regularly_sumstats")

results_df_other <- map_dfr(trait_list_other, my_function)

results_df_other <- results_df_other %>%
  mutate(P_adjusted = p.adjust(p, method = "BH", n = 11)) 

write_delim(results_df_other, file = "Compare_female_rg_to_male_rg_for_other_traits.txt")




trait_list_BMI_sex_specific <- c("BMI_female_sumstats", "BMI_male_sumstats")

results_df_BMI_sex_specific <- map_dfr(trait_list_BMI_sex_specific, my_function)

results_df_BMI_sex_specific <- results_df_BMI_sex_specific %>%
  mutate(P_adjusted = p.adjust(p, method = "BH", n = 2))


write_delim(results_df_BMI_sex_specific, file = "Compare_female_rg_to_male_rg_for_BMI_sex_specific.txt")















