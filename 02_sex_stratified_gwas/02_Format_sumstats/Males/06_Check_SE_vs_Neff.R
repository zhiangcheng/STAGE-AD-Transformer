# Plot mean of SE vs effective sample size of each cohort - expect SE to reflect Neff

directory = "/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/Males/"

GWAS_directory = "/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/010_Format_sumstats/"

ADNI_male = paste0(GWAS_directory, "Males/adni_ad_male_sumstats_formatted_formetaanalysis.txt", sep = "")

MSBB_male = paste0(GWAS_directory, "Males/msbb_ad_male_sumstats_formatted_formetaanalysis.txt", sep = "")

ROSMAP_male = paste0(GWAS_directory, "Males/rosmap_ad_male_sumstats_formatted_formetaanalysis.txt", sep = "")


## ---- LoadPackages
library(tidyverse)
## ----

## ---- LoadData
load_dataframe <- function(df_name) {
  df <- data.table::fread(df_name, header = T, stringsAsFactors = F)
  
  return(df)
}

file_paths <- c(
  ADNI_male = ADNI_male,
  MSBB_male = MSBB_male,
  ROSMAP_male = ROSMAP_male
)

data_frames <- lapply(file_paths, load_dataframe)
# Access individual data frames
# data_frames[["ADNI_male"]]
## ----

## ---- Dataframe
# Calculate Neff and mean of beta SE for each cohort using just SNPs with freq of 0.45 - 0.55

process_dataset <- function(dataset, cohort_name, N, P) {
  dataset %>%
    summarise(Cohort = cohort_name, 
              Mean_SE_0_0.2 = mean(SE[FREQA1 > 0 & FREQA1 < 0.2], na.rm = TRUE),
              Mean_SE_0.2_0.45 = mean(SE[FREQA1 > 0.2 & FREQA1 < 0.45], na.rm = TRUE),
              Mean_SE_0.45_0.55 = mean(SE[FREQA1 > 0.45 & FREQA1 < 0.55], na.rm = TRUE),
              Mean_SE_all = mean(SE, na.rm = TRUE)) %>%
    mutate(Neff = 4 * N * P * (1 - P)) %>%
    pivot_longer(cols = starts_with("Mean_SE"), 
                 names_to = "FREQA1_Bracket", 
                 values_to = "Mean_SE") %>%
    mutate(FREQA1_Bracket = case_when(
      FREQA1_Bracket == "Mean_SE_0_0.2" ~ "0-0.2",
      FREQA1_Bracket == "Mean_SE_0.2_0.45" ~ "0.2-0.45",
      FREQA1_Bracket == "Mean_SE_0.45_0.55" ~ "0.45-0.55",
      FREQA1_Bracket == "Mean_SE_all" ~ "All"
    ))
  }

# Define cohort name, N and P for each dataset/cohort
cohort_details <- list(
  ADNI = list(dataset = data_frames[["ADNI_male"]], N = 314, P = 0.2803),
  MSBB = list(dataset = data_frames[["MSBB_male"]], N = 226, P = 0.5044),
  ROSMAP = list(dataset = data_frames[["ROSMAP_male"]], N = 747, P = 0.4485)
)


result <- lapply(names(cohort_details), function(cohort_name) {
  cohort <- cohort_details[[cohort_name]] # Extract details for the current cohort
  process_dataset(
    dataset = cohort$dataset,
    cohort_name = cohort_name,
    N = cohort$N,
    P = cohort$P
  )
})


# Combine the results into a single dataframe
summary_df <- do.call(rbind, result)
## ----

## ---- Plot
plot <- ggplot(summary_df, aes(x = Neff, y = Mean_SE)) +
  geom_point() +
  geom_text(aes(label = Cohort), nudge_y = 0.002) + 
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ FREQA1_Bracket, scales = "free") +
  scale_x_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, 200),
                     labels = seq(0, 1000, 200)) +
  ggtitle("Mean SE vs Neff for SNPs by FREQA1 Bracket") +
  theme_classic()

plot

outfile = paste0(directory, "SE_vs_Neff_male_cohorts.pdf")
ggsave(plot, file = outfile, width = 30, height = 30, unit = "cm")
## ----
