# This script tests if female and male betas are different
# Does regression of female-female, male-male and female-male betas for SNPs known to be associated with depression in each cohort
# Then meta-analyses each of these regression betas across all the cohorts
# Across sexes there are 30 pairwise linear regression to meta-analyse
# Within the same sex as linear regression is directional can only meta-analyse 15 linear regressions (A vs B and B vs A will have a different slope but they are not independent so can't go in the same meta-analysis)
# There are 32,768 different sets of 15 linear regressions that can be meta-analysed - so do all of these and calculate the mean and SD then compare to the male vs female comparison.

directory = "/path/09_Effect_sizes_plots/Adams_SNPs/"

Adams_clump = "/path/09_Effect_sizes_plots/Adams_SNPs/clumped_Adams_SNPs_all_results.txt"

GWAS_directory = "/path/01_Format_sumstats/"

AGDS_female = paste0(GWAS_directory, "Females/AGDS_MDD_female_sumstats_formatted_formetaanalysis.txt", sep = "")

AllOfUs_female = paste0(GWAS_directory, "Females/AllOfUs_MDD_female_sumstats_formatted_formetaanalysis.txt", sep = "")

Bionic_female = paste0(GWAS_directory, "Females/Bionic_MDD_female_sumstats_formatted_formetaanalysis.txt", sep = "")

Blokland_female = paste0(GWAS_directory, "Females/Blokland21_MDD_female_sumstats_formatted_formetaanalysis.txt")

GLAD_female = paste0(GWAS_directory, "Females/GLAD_MDD_female_sumstats_formatted_formetaanalysis.txt", sep = "")

UKB_female = paste0(GWAS_directory, "Females/UKB_MDD_female_sumstats_formatted_formetaanalysis.txt", sep = "")


AGDS_male = paste0(GWAS_directory, "Males/AGDS_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

AllOfUs_male = paste0(GWAS_directory, "Males/AllOfUs_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

Bionic_male = paste0(GWAS_directory, "Males/Bionic_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

Blokland_male = paste0(GWAS_directory, "Males/Blokland21_MDD_male_sumstats_formatted_formetaanalysis.txt")

GLAD_male = paste0(GWAS_directory, "Males/GLAD_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

UKB_male = paste0(GWAS_directory, "Males/UKB_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")


GxS = "/path/03_Metal/GxS/fulldc/Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID.txt"


MA_female = "/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

MA_male = "/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"

## ---- LoadPackages
library(performance)
library(DHARMa)
library(ggstance)   # for horizontal error bars
library(emmeans)
library(ggpubr)     # add equation to plot
library(metafor)    # for meta-regression
library(gridExtra)  # for grid layout of plots
library(broom)      # for tidy outputs
library(scales)     # remove scientific notation from ggplot axis labels
library(tidyverse)
## ----


## ---- LoadData
load_dataframe <- function(df_name) {
  df <- read.table(df_name, header = T, stringsAsFactors = F)
  
  return(df)
}

file_paths <- c(
  AGDS_female = AGDS_female,
  AllOfUs_female = AllOfUs_female,
  Bionic_female = Bionic_female,
  Blokland_female = Blokland_female,
  GLAD_female = GLAD_female,
  UKB_female = UKB_female,
  AGDS_male = AGDS_male,
  AllOfUs_male = AllOfUs_male,
  Bionic_male = Bionic_male,
  Blokland_male = Blokland_male,
  GLAD_male = GLAD_male,
  UKB_male = UKB_male
)

data_frames <- lapply(file_paths, load_dataframe)

# Access individual data frames
# data_frames[["AGDS_female"]]

Adams_clump_df <- read.table(Adams_clump, header = T, stringsAsFactors = F)

GxS_df <- read.table(GxS, header = T, stringsAsFactors = F)

MA_female_df <- read.table(MA_female, header = T)

MA_male_df <- read.table(MA_male, header = T)
## ----


## ---- AlignSNPs
# Pull out independent lead SNPs genome-wide sig from Adams
Adams_clump_df_format <- Adams_clump_df %>%
  filter(P < 5e-08)

Adams_sig_SNPs <- Adams_clump_df_format %>%
  pull(SNP)

# Choose tested allele for each SNP to use in each cohort - based on minor allele in GxS meta-analysis sumstats
consistent_alleles <- GxS_df %>%
  filter(rsID_build37 %in% Adams_sig_SNPs) %>%
  mutate(Consistent_A1 = case_when(
    Freq1 < 0.5 ~ Allele1,
    Freq1 > 0.5 ~ Allele2),
    Consistent_A2 =  case_when(
      Freq1 < 0.5 ~ Allele2,
      Freq1 > 0.5 ~ Allele1)) %>%
  rename(MARKER_build37 = MarkerName) %>%
  select(MARKER_build37, Consistent_A1, Consistent_A2) %>%
  mutate(Consistent_A1 = toupper(Consistent_A1),
         Consistent_A2 = toupper(Consistent_A2))

# Now restrict each sumstats to Adams lead independent genome-wide significant SNPs and make sure tested allele consistent for all sumstats
process_and_align_dataset <- function(dataset, cohort_name, sex) {
  dataset <- dataset %>%
    inner_join(consistent_alleles, by = "MARKER_build37") %>%
    mutate(
      Flip = (A1 != Consistent_A1),
      BETA = ifelse(Flip, -BETA, BETA),
      A1 = ifelse(Flip, Consistent_A1, A1),
      A2 = ifelse(Flip, Consistent_A2, A2)
    ) %>%
    select(-Flip, -Consistent_A1, -Consistent_A2) %>%
    rename_with(~paste0(., "_", cohort_name, "_", sex), -MARKER_build37)
  
  return(dataset)
}

female_datasets <- list(
  AGDS = data_frames[["AGDS_female"]],
  AllOfUs = data_frames[["AllOfUs_female"]],
  Bionic = data_frames[["Bionic_female"]],
  Blokland = data_frames[["Blokland_female"]],
  GLAD = data_frames[["GLAD_female"]],
  UKB = data_frames[["UKB_female"]]
)

male_datasets <- list(
  AGDS = data_frames[["AGDS_male"]],
  AllOfUs = data_frames[["AllOfUs_male"]],
  Bionic = data_frames[["Bionic_male"]],
  Blokland = data_frames[["Blokland_male"]],
  GLAD = data_frames[["GLAD_male"]],
  UKB = data_frames[["UKB_male"]]
)

all_datasets <- list(
  female = female_datasets,
  male = male_datasets
)


# Process all datasets and add suffixes
aligned_datasets <- lapply(names(all_datasets), function(sex) {
  datasets <- all_datasets[[sex]]
  lapply(names(datasets), function(cohort_name) {
    process_and_align_dataset(datasets[[cohort_name]], cohort_name, sex)
  })
})

# Flatten the nested list into a single list of data frames
aligned_datasets_flat <- unlist(aligned_datasets, recursive = FALSE)

# Combine all datasets into one dataframe by joining on MARKER_build37
processed_dataframes <- Reduce(function(x, y) full_join(x, y, by = "MARKER_build37"), aligned_datasets_flat)


outfile <- paste0(directory, "/processed_dataframes.RData", sep = "")
save(processed_dataframes, file = outfile)
## ----


## ---- CorrelationFunctions
# standarized beta and SE in SD unit
# z represent z statistics, p presents allele frequency, n represents sample size.\
calcu_std_b_se<-function(z,p,n){
  std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
  std_se=1/sqrt(2*p*(1-p)*(n+z^2))
  res<-data.frame(std_b_hat,std_se);
  return(res)
}
## ----


## ---- CorrelationPlotRegressionMalevsFemale
generate_plots_and_save_regressions <- function(processed_dataframes) {
  # Define sample sizes for each cohort and sex
  sample_sizes <- list(
    AGDS = list(female = 17553, male = 9775),
    AllOfUs = list(female = 73454, male = 46972),
    Bionic = list(female = 37542, male = 24445),
    Blokland = list(female = 41937, male = 28993),
    GLAD = list(female = 20101, male = 7681),
    UKB = list(female = 99405, male = 79124)
  )
  
  cohort_names <- names(sample_sizes)
  
  # Generate all combinations of female and male cohorts
  combinations <- expand.grid(female_cohort = cohort_names, male_cohort = cohort_names)
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Create a data frame to store results
  regression_results <- data.frame(
    female_cohort = character(),
    male_cohort = character(),
    regression_slope = numeric(),
    regression_slope_se = numeric(),
    regression_slope_p = numeric(),
    regression_intercept = numeric(),
    regression_intercept_se = numeric(),
    regression_intercept_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each combination
  for (i in seq_len(nrow(combinations))) {
    female_cohort <- combinations$female_cohort[i]
    male_cohort <- combinations$male_cohort[i]
    
    # Extract relevant columns for the current combination
    data <- processed_dataframes %>%
      select(
        MARKER_build37,
        starts_with(paste0("BETA_", female_cohort, "_female")),
        starts_with(paste0("SE_", female_cohort, "_female")),
        starts_with(paste0("BETA_", male_cohort, "_male")),
        starts_with(paste0("SE_", male_cohort, "_male")),
        starts_with(paste0("FREQA1_", female_cohort, "_female")),
        starts_with(paste0("FREQA1_", male_cohort, "_male"))
      ) %>%
      rename_with(~ sub(paste0("_", female_cohort, "_female"), "_female", .), starts_with("BETA")) %>%
      rename_with(~ sub(paste0("_", male_cohort, "_male"), "_male", .), starts_with("BETA")) %>%
      rename_with(~ sub(paste0("_", female_cohort, "_female"), "_female", .), starts_with("SE")) %>%
      rename_with(~ sub(paste0("_", male_cohort, "_male"), "_male", .), starts_with("SE")) %>%
      rename_with(~ sub(paste0("_", female_cohort, "_female"), "_female", .), starts_with("FREQA1")) %>%
      rename_with(~ sub(paste0("_", male_cohort, "_male"), "_male", .), starts_with("FREQA1"))
    
    # Standardize betas and SEs for female and male cohorts
    female_std <- calcu_std_b_se(
      z = data$BETA_female / data$SE_female,
      p = data$FREQA1_female,
      n = sample_sizes[[female_cohort]]$female
    )
    male_std <- calcu_std_b_se(
      z = data$BETA_male / data$SE_male,
      p = data$FREQA1_male,
      n = sample_sizes[[male_cohort]]$male
    )
    
    # Create a new dataframe with standardized betas and standard errors
    standardized_data <- data.frame(
      MARKER_build37 = data$MARKER_build37,
      female_std_b = female_std$std_b_hat,
      female_std_se = female_std$std_se,
      male_std_b = male_std$std_b_hat,
      male_std_se = male_std$std_se
    )
    
    
    # Run linear regression of standardized male beta vs. standardized female beta
    lm_model <- lm(male_std_b ~ female_std_b, data = standardized_data)
    
    model_summary <- summary(lm_model)
    
    # Extract regression coefficients
    regression_slope <- model_summary$coefficients["female_std_b", "Estimate"]
    regression_slope_se <- model_summary$coefficients["female_std_b", "Std. Error"]
    regression_slope_p <- model_summary$coefficients["female_std_b", "Pr(>|t|)"]
    regression_intercept <- model_summary$coefficients["(Intercept)", "Estimate"]
    regression_intercept_se <- model_summary$coefficients["(Intercept)", "Std. Error"]
    regression_intercept_p <- model_summary$coefficients["(Intercept)", "Pr(>|t|)"]
    
    
    # Add results to the regression_results data frame
    regression_results <- rbind(
      regression_results,
      data.frame(
        female_cohort = female_cohort,
        male_cohort = male_cohort,
        regression_slope = regression_slope,
        regression_slope_se = regression_slope_se,
        regression_slope_p = regression_slope_p,
        regression_intercept = regression_intercept,
        regression_intercept_se = regression_intercept_se,
        regression_intercept_p = regression_intercept_p
      )
    )
    
    # Calculate axis limits based on standardized values
    x_min <- min(standardized_data$female_std_b - standardized_data$female_std_se, na.rm = TRUE)
    x_max <- max(standardized_data$female_std_b + standardized_data$female_std_se, na.rm = TRUE)
    y_min <- min(standardized_data$male_std_b - standardized_data$male_std_se, na.rm = TRUE)
    y_max <- max(standardized_data$male_std_b + standardized_data$male_std_se, na.rm = TRUE)
    axis_limit <- max(abs(c(x_min, x_max, y_min, y_max)))
    label_x <- -0.9 * axis_limit
    label_y <- 0.9 * axis_limit
    
    # Generate the plot using standardized betas
    plot <- ggplot(standardized_data, aes(x = female_std_b, y = male_std_b)) +
      geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
      geom_linerange(aes(ymin = male_std_b - male_std_se, ymax = male_std_b + male_std_se), size = 0.4) +
      geom_linerangeh(aes(xmin = female_std_b - female_std_se, xmax = female_std_b + female_std_se), size = 0.4) +
      geom_point(alpha = 0.75) +
      geom_abline(colour = "grey40") +
      geom_smooth(method = "lm", se = FALSE, color = "grey80") +
      scale_x_continuous(paste("Standardized Female Beta (", female_cohort, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      scale_y_continuous(paste("Standardized Male Beta (", male_cohort, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      annotate("text", x = label_x, y = label_y,
               label = paste0(
                 "Intercept = ", round(regression_intercept, 3), " ± ", round(regression_intercept_se, 3),
                 ", P = ", format.pval(regression_intercept_p, digits = 3, scientific = TRUE), "\n",
                 "Beta = ", round(regression_slope, 3), " ± ", round(regression_slope_se, 3),
                 ", P = ", format.pval(regression_slope_p, digits = 3, scientific = TRUE)
               ),
               hjust = 0) +
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(female_cohort, male_cohort, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, regressions = regression_results))
}


# Generate all plots
all_plots_regressions_male_v_female <- generate_plots_and_save_regressions(processed_dataframes)

# Example: Display one plot
print(all_plots_regressions_male_v_female[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display regression results
all_plots_regressions_male_v_female$regressions


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_regressions_male_v_female$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/male_vs_female_betas_linear_regression_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----


## ---- MetaAnalysisMaleVsFemale
# Meta-analyse (inverse variance weighted) the correlations and betas for the M-F regression in each of the 6 cohorts
# Uses random effects (heterogeneity between cohorts)

regression_male_v_female_df <- all_plots_regressions_male_v_female$regressions %>%
  as_tibble() %>%
  mutate(cohort_label = paste(female_cohort, "vs", male_cohort))

outfile = paste0(directory, "/male_vs_female_linear_regression_results_all_cohort_combos.RData")
save(regression_male_v_female_df, file = outfile)

## Across cohorts (all comparisons apart from diagonal) = 30 comparisons

regression_male_v_female_across_df <- regression_male_v_female_df %>%
  filter(female_cohort != male_cohort)

# Meta-analysis of male beta vs female beta across cohorts regression slope (beta value from regression)
meta_slope_MF_across <- rma(yi = regression_slope,
                            sei = regression_slope_se,
                            weighted = T,
                            data = regression_male_v_female_across_df)


summary(meta_slope_MF_across)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_slope_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_slope_MF_across, slab = regression_male_v_female_across_df$cohort_label, xlab = "Slope of male beta vs female beta")
dev.off()


# Meta-analysis of male beta vs female beta across cohorts regression intercept
meta_intcpt_MF_across <- rma(yi = regression_intercept,
                             sei = regression_intercept_se,
                             weighted = T,
                             data = regression_male_v_female_across_df)




# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_intercept_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_intcpt_MF_across, slab = regression_male_v_female_across_df$cohort_label, xlab = "Intercept of male beta vs female beta")
dev.off()


# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_linear_regression_summary.txt")
sink(outfile)

## Slope
# Write the summary of meta_slope
cat("\n\nSummary of meta-analysis of slope (beta value from regression of male beta vs female beta across cohorts):\n")
print(summary(meta_slope_MF_across))

# estimate significantly different to 1?
estimate <- summary(meta_slope_MF_across)[["beta"]]
SE <- summary(meta_slope_MF_across)[["se"]]
hypothesized_value <- 1

Z <- (estimate - hypothesized_value) / SE
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = TRUE)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 1):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 1 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 1 (p >= 0.05).\n")
}


## Intercept
# Write the summary of meta_intcpt
cat("\n\nSummary of meta-analysis of intercept (intercept value from regression of male beta vs female beta across cohorts):\n")
print(summary(meta_intcpt_MF_across))

# estimate significantly different to 0?
estimate <- summary(meta_intcpt_MF_across)[["beta"]]
SE <- summary(meta_intcpt_MF_across)[["se"]]
hypothesized_value <- 0

Z <- (estimate - hypothesized_value) / SE
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = F)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 0):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 0 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 0 (p >= 0.05).\n")
}


# Stop capturing output
sink()




## Within cohorts (all diagonal comparisons)

regression_male_v_female_within_df <- regression_male_v_female_df %>%
  filter(female_cohort == male_cohort)

# Meta-analysis of male beta vs female beta within same cohort regression slope (beta value from regression)
meta_slope_MF_within <- rma(yi = regression_slope,
                            sei = regression_slope_se,
                            weighted = T,
                            data = regression_male_v_female_within_df)


summary(meta_slope_MF_within)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_slope_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_slope_MF_within, slab = regression_male_v_female_within_df$cohort_label, xlab = "Slope of male beta vs female beta")
dev.off()


# Meta-analysis of male beta vs female beta within the same cohort regression intercept
meta_intcpt_MF_within <- rma(yi = regression_intercept,
                             sei = regression_intercept_se,
                             weighted = T,
                             data = regression_male_v_female_within_df)


summary(meta_intcpt_MF_within)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_intercept_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_intcpt_MF_within, slab = regression_male_v_female_within_df$cohort_label, xlab = "Intercept of male beta vs female beta")
dev.off()


# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_linear_regression_summary.txt")
sink(outfile)

## Slope
# Write the summary of meta_slope
cat("\n\nSummary of meta-analysis of slope (beta value from regression of male beta vs female beta within the same cohort):\n")
print(summary(meta_slope_MF_within))

# estimate significantly different to 1?
estimate <- summary(meta_slope_MF_within)[["beta"]]
SE <- summary(meta_slope_MF_within)[["se"]]
hypothesized_value <- 1

Z <- (estimate - hypothesized_value) / SE
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = TRUE)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 1):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 1 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 1 (p >= 0.05).\n")
}


## Intercept
# Write the summary of meta_intcpt
cat("\n\nSummary of meta-analysis of intercept (intercept value from regression of male beta vs female beta within the same cohort):\n")
print(summary(meta_intcpt_MF_within))

# estimate significantly different to 0?
estimate <- summary(meta_intcpt_MF_within)[["beta"]]
SE <- summary(meta_intcpt_MF_within)[["se"]]
hypothesized_value <- 0

Z <- (estimate - hypothesized_value) / SE
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = F)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 0):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 0 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 0 (p >= 0.05).\n")
}


# Stop capturing output
sink()
## ----


## ---- CorrelationPlotRegressionMalevsMale
generate_plots_and_save_male_male_regression <- function(processed_dataframes) {
  # Define sample sizes for each cohort and sex
  sample_sizes <- list(
    AGDS = list(male = 9775),
    AllOfUs = list(male = 46972),
    Bionic = list(male = 24445),
    Blokland = list(male = 28993),
    GLAD = list(male = 7681),
    UKB = list(male = 79124)
  )
  
  cohort_names <- names(sample_sizes)
  
  # Generate all combinations of male cohorts (for male vs male comparison)
  combinations <- expand.grid(male_cohort1 = cohort_names, male_cohort2 = cohort_names)
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Create a data frame to store correlation results for male vs male
  regression_results <- data.frame(
    male_cohort1 = character(),
    male_cohort2 = character(),
    regression_slope = numeric(),
    regression_slope_se = numeric(),
    regression_slope_p = numeric(),
    regression_intercept = numeric(),
    regression_intercept_se = numeric(),
    regression_intercept_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each combination
  for (i in seq_len(nrow(combinations))) {
    male_cohort1 <- combinations$male_cohort1[i]
    male_cohort2 <- combinations$male_cohort2[i]
    
    # Extract relevant columns for the current combination
    if (male_cohort1 == male_cohort2) {
      data <- processed_dataframes %>%
        select(
          MARKER_build37,
          starts_with(paste0("BETA_", male_cohort1, "_male")),
          starts_with(paste0("SE_", male_cohort1, "_male")),
          starts_with(paste0("FREQA1_", male_cohort1, "_male"))
        ) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("FREQA1")) %>%
        mutate(
          BETA_male2 = BETA_male1,  # Copy BETA_male1 to BETA_male2
          SE_male2 = SE_male1,      # Copy SE_male1 to SE_male2
          FREQA1_male2 = FREQA1_male1  # Copy FREQA1_male1 to FREQA1_male2
        )
    } else {
      data <- processed_dataframes %>%
        select(
          MARKER_build37,
          starts_with(paste0("BETA_", male_cohort1, "_male")),
          starts_with(paste0("SE_", male_cohort1, "_male")),
          starts_with(paste0("BETA_", male_cohort2, "_male")),
          starts_with(paste0("SE_", male_cohort2, "_male")),
          starts_with(paste0("FREQA1_", male_cohort1, "_male")),
          starts_with(paste0("FREQA1_", male_cohort2, "_male"))
        ) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", male_cohort2, "_male"), "_male2", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", male_cohort2, "_male"), "_male2", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", male_cohort1, "_male"), "_male1", .), starts_with("FREQA1")) %>%
        rename_with(~ sub(paste0("_", male_cohort2, "_male"), "_male2", .), starts_with("FREQA1"))
    }
    
    
    # Standardize betas and SEs for both male cohorts
    male1_std <- calcu_std_b_se(
      z = data$BETA_male1 / data$SE_male1,
      p = data$FREQA1_male1,
      n = sample_sizes[[male_cohort1]]$male
    )
    male2_std <- calcu_std_b_se(
      z = data$BETA_male2 / data$SE_male2,
      p = data$FREQA1_male2,
      n = sample_sizes[[male_cohort2]]$male
    )
    
    # Create a new dataframe with standardized betas and standard errors
    standardized_data <- data.frame(
      MARKER_build37 = data$MARKER_build37,
      male1_std_b = male1_std$std_b_hat,
      male1_std_se = male1_std$std_se,
      male2_std_b = male2_std$std_b_hat,
      male2_std_se = male2_std$std_se
    )
    
    
    # Run linear regression of standardized male1 beta vs. standardized male2 beta
    lm_model <- lm(male2_std_b ~ male1_std_b, data = standardized_data)
    
    model_summary <- summary(lm_model)
    
    # Extract regression coefficients
    regression_slope <- model_summary$coefficients["male1_std_b", "Estimate"]
    regression_slope_se <- model_summary$coefficients["male1_std_b", "Std. Error"]
    regression_slope_p <- model_summary$coefficients["male1_std_b", "Pr(>|t|)"]
    regression_intercept <- model_summary$coefficients["(Intercept)", "Estimate"]
    regression_intercept_se <- model_summary$coefficients["(Intercept)", "Std. Error"]
    regression_intercept_p <- model_summary$coefficients["(Intercept)", "Pr(>|t|)"]
    
    # Add results to the regression_results data frame
    regression_results <- rbind(
      regression_results,
      data.frame(
        male_cohort1 = male_cohort1,
        male_cohort2 = male_cohort2,
        regression_slope = regression_slope,
        regression_slope_se = regression_slope_se,
        regression_slope_p = regression_slope_p,
        regression_intercept = regression_intercept,
        regression_intercept_se = regression_intercept_se,
        regression_intercept_p = regression_intercept_p
      )
    )
    
    # Calculate axis limits based on standardized values
    x_min <- min(standardized_data$male1_std_b - standardized_data$male1_std_se, na.rm = TRUE)
    x_max <- max(standardized_data$male1_std_b + standardized_data$male1_std_se, na.rm = TRUE)
    y_min <- min(standardized_data$male2_std_b - standardized_data$male2_std_se, na.rm = TRUE)
    y_max <- max(standardized_data$male2_std_b + standardized_data$male2_std_se, na.rm = TRUE)
    axis_limit <- max(abs(c(x_min, x_max, y_min, y_max)))
    label_x <- -0.9 * axis_limit
    label_y <- 0.9 * axis_limit
    
    # Generate the plot using standardized betas
    plot <- ggplot(standardized_data, aes(x = male1_std_b, y = male2_std_b)) +
      geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
      geom_linerange(aes(ymin = male2_std_b - male2_std_se, ymax = male2_std_b + male2_std_se), size = 0.4) +
      geom_linerangeh(aes(xmin = male1_std_b - male1_std_se, xmax = male1_std_b + male1_std_se), size = 0.4) +
      geom_point(alpha = 0.75) +
      geom_abline(colour = "grey40") +
      geom_smooth(method = "lm", se = FALSE, color = "grey80") +
      scale_x_continuous(paste("Standardized Male Beta (", male_cohort1, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      scale_y_continuous(paste("Standardized Male Beta (", male_cohort2, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      annotate("text", x = label_x, y = label_y,
               label = paste0(
                 "Intercept = ", round(regression_intercept, 3), " ± ", round(regression_intercept_se, 3),
                 ", P = ", format.pval(regression_intercept_p, digits = 3, scientific = TRUE), "\n",
                 "Beta = ", round(regression_slope, 3), " ± ", round(regression_slope_se, 3),
                 ", P = ", format.pval(regression_slope_p, digits = 3, scientific = TRUE)
               ),
               hjust = 0) +
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(male_cohort1, male_cohort2, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, regression = regression_results))
}


# Generate all plots
all_plots_regressions_male_v_male <- generate_plots_and_save_male_male_regression(processed_dataframes)

# Example: Display one plot
print(all_plots_regressions_male_v_male[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display correlations
all_plots_regressions_male_v_male$regression


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_regressions_male_v_male$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/male_vs_male_betas_linear_regression_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----


## ---- MetaAnalysisMaleVsMale

regression_male_v_male_df <- all_plots_regressions_male_v_male$regression %>%
  as_tibble() %>%
  mutate(cohort_label = paste(male_cohort1, "vs", male_cohort2))


# Create unique unordered pairs (15 independent pairs)
cohorts <- unique(c(regression_male_v_male_df$male_cohort1, regression_male_v_male_df$male_cohort2))
independent_pairs <- combn(cohorts, 2, simplify = FALSE)  # 15 unique unordered pairs

# Create both directions for each pair (A vs B and B vs A)
pairs_directions <- lapply(independent_pairs, function(pair) {
  list(paste(pair[1], pair[2], sep = " vs "), paste(pair[2], pair[1], sep = " vs "))
})

# Generate all possible sets of 15 pairs by selecting one direction per pair
# Each of the 15 pairs has 2 possible directions, so there are 2^15 = 32,768 sets
direction_combinations <- expand.grid(rep(list(c(1, 2)), 15))  # 32,768 rows, each a unique combination of directions (1st or 2nd option for each combo)

# Meta-analyse SLOPE
meta_results_slope <- list()

for (i in seq_len(nrow(direction_combinations))) {
  selected_pairs <- sapply(seq_along(pairs_directions), function(j) {
    pairs_directions[[j]][[direction_combinations[i, j]]]  # Pick A->B or B->A for each pair
  })
  
  # Filter the regression dataset to include only selected pairs
  selected_data <- regression_male_v_male_df %>%
    filter(cohort_label %in% selected_pairs)
  
  # Run the meta-analysis
  meta_results_slope[[paste("Meta_analysis_set", i, sep = "_")]] <- rma(
    yi = selected_data$regression_slope,
    sei = selected_data$regression_slope_se,
    weighted = TRUE,
    data = selected_data
  )
}

# Meta-analyse INTERCEPT
meta_results_intcpt <- list()

for (i in seq_len(nrow(direction_combinations))) {
  selected_pairs <- sapply(seq_along(pairs_directions), function(j) {
    pairs_directions[[j]][[direction_combinations[i, j]]]  # Pick A->B or B->A for each pair
  })
  
  # Filter the regression dataset to include only selected pairs
  selected_data <- regression_male_v_male_df %>%
    filter(cohort_label %in% selected_pairs)
  
  # Run the meta-analysis
  meta_results_intcpt[[paste("Meta_analysis_set", i, sep = "_")]] <- rma(
    yi = selected_data$regression_intercept,
    sei = selected_data$regression_intercept_se,
    weighted = TRUE,
    data = selected_data
  )
}

# View results
summary(meta_results_slope[[1]])  # View first meta-analysis result SLOPE
summary(meta_results_intcpt[[1]])  # View first meta-analysis result INTERCEPT

# Save results
# SLOPE
meta_analyses_MM_slope_df <- do.call(rbind, lapply(seq_along(meta_results_slope), function(i) {
  meta <- meta_results_slope[[i]]
  
  # Extract relevant values from the model summary
  data.frame(
    meta_set = paste("Meta_analysis_set", i, sep = "_"),
    tau2 = meta$tau2,   # Between-study variance
    tau = sqrt(meta$tau2),
    I2 = meta$I2,       # Percentage of variability due to heterogeneity
    H2 = meta$H2,       # Total variability / sampling variability
    Q = meta$QE,        # Cochran’s Q statistic
    Q_pval = meta$QEp,  # p-value for heterogeneity test
    estimate = meta$b,  # Meta-analysis estimated effect size
    se = meta$se,       # Standard error of the effect size
    sd = meta$se * sqrt(15),  # Standard deviation (SD = SE * sqrt(no. studies included = 15 in each))
    zval = meta$zval,   # Z-statistic
    pval = meta$pval,   # p-value for meta-analysis effect
    ci_lb = meta$ci.lb, # Lower bound of 95% CI
    ci_ub = meta$ci.ub  # Upper bound of 95% CI
  )
}))


outfile <- paste0(directory, "male_vs_male_meta_analyses_slope_all_32768_sets_of_15_linear_regressions_results.RData", sep = "")
save(meta_analyses_MM_slope_df, file = outfile)


# INTERCEPT
meta_analyses_MM_intcpt_df <- do.call(rbind, lapply(seq_along(meta_results_intcpt), function(i) {
  meta <- meta_results_intcpt[[i]]
  
  # Extract relevant values from the model summary
  data.frame(
    meta_set = paste("Meta_analysis_set", i, sep = "_"),
    tau2 = meta$tau2,   # Between-study variance
    tau = sqrt(meta$tau2),
    I2 = meta$I2,       # Percentage of variability due to heterogeneity
    H2 = meta$H2,       # Total variability / sampling variability
    Q = meta$QE,        # Cochran’s Q statistic
    Q_pval = meta$QEp,  # p-value for heterogeneity test
    estimate = meta$b,  # Meta-analysis estimated effect size
    se = meta$se,       # Standard error of the effect size
    sd = meta$se * sqrt(15),  # Standard deviation (SD = SE * sqrt(no. studies included = 15 in each))
    zval = meta$zval,   # Z-statistic
    pval = meta$pval,   # p-value for meta-analysis effect
    ci_lb = meta$ci.lb, # Lower bound of 95% CI
    ci_ub = meta$ci.ub  # Upper bound of 95% CI
  )
}))


outfile <- paste0(directory, "male_vs_male_meta_analyses_intcpt_all_32768_sets_of_15_linear_regressions_results.RData", sep = "")
save(meta_analyses_MM_intcpt_df, file = outfile)

# Calculate mean of all meta-analyses estimates, and sd of this mean estimate
MM_slope_means <- meta_analyses_MM_slope_df %>%
  summarise(mean_estimate = mean(estimate),
            mean_estimate_sd = sd(estimate),
            mean_estimate_se = mean_estimate_sd / sqrt(n()))

hist(meta_analyses_MM_slope_df$estimate)


MM_intcpt_means <- meta_analyses_MM_intcpt_df %>%
  summarise(mean_estimate = mean(estimate),
            mean_estimate_sd = sd(estimate),
            mean_estimate_se = mean_estimate_sd / sqrt(n()))

hist(meta_analyses_MM_intcpt_df$estimate)



# Define output file
outfile <- paste0(directory, "/meta_male_vs_male_across_cohorts_linear_regression_summary.txt")
sink(outfile)

## Slope
# Write the mean and sd of slope distributions
cat("\n\nResults from slope (beta value from regression of male beta vs male beta across cohorts) were meta-analysed in sets of 15 (32,768 meta-analyses):\n")

MM_slope_mean <- MM_slope_means %>% pull(mean_estimate)
MM_slope_sd <- MM_slope_means %>% pull(mean_estimate_sd)
MM_slope_se <- MM_slope_means %>% pull(mean_estimate_se)

cat("Mean slope from distribution of meta-anlyses", MM_slope_mean, "\n")
cat("SD of mean slope from distribution of meta-anlyses", MM_slope_sd, "\n")
cat("SE of mean slope from distribution of meta-anlyses", MM_slope_se, "\n")


## Intercept
cat("\n\nResults from intercept (beta value from regression of male beta vs male beta across cohorts) were meta-analysed in sets of 15 (32,768 meta-analyses):\n")

MM_intercept_mean <- MM_intcpt_means %>% pull(mean_estimate)
MM_intercept_sd <- MM_intcpt_means %>% pull(mean_estimate_sd)
MM_intercept_se <- MM_intcpt_means %>% pull(mean_estimate_se)

cat("Mean intercept from distribution of meta-anlyses", MM_intercept_mean, "\n")
cat("SD of mean intercept from distribution of meta-anlyses", MM_intercept_sd, "\n")
cat("SE of mean intercept from distribution of meta-anlyses", MM_intercept_se, "\n")

# Stop capturing output
sink()
## ----






## ---- CompareMFvsMM

## SLOPE

# Does the observed M vs. F slope differ significantly from the distribution of M vs. M slopes
# Z = beta(male vs female meta- analysis) - mean (of all male vs male meta-analysis estimates) / sd (of all male vs male meta-analysis estimates)
# Then calculate p-value from Z-score
# If p < 0.05 = male vs female slope is sig different to baseline male vs male slope


Z_MF_MM_slope <- (meta_slope_MF_across$b - MM_slope_means$mean_estimate) / MM_slope_means$mean_estimate_sd
p_MF_MM_slope <- 2 * pnorm(abs(Z_MF_MM_slope), mean = 0, sd = 1, lower.tail = F)



# Visualise
# Histogram of M vs M estimates and overlay M vs F slope
plot_distr_slope <- ggplot(meta_analyses_MM_slope_df, aes(x = estimate)) +
  geom_histogram(binwidth = 0.005, fill = "grey", colour = "black") +
  geom_vline(xintercept = meta_slope_MF_across$b[1,1], colour = "red") +
  geom_text(label = "Meta-analysis \nestimate of \nMale vs Female \nlinear regression \nslope",
            x = 0.39, y = 1500, colour = "red", hjust = 0, vjust = 1) +
  scale_x_continuous("Meta-analysis estimates of Male vs Male linear regression slope") +
  scale_y_continuous("Count") +
  theme_classic()

outfile = paste0(directory, "/male_vs_male_slopes_distribution_overlaid_MF_slope.png")
ggsave(plot_distr_slope, file = outfile, width = 30, height = 30, unit = "cm")




## INTERCEPT

# Does the observed M vs. F intcpt differ significantly from the distribution of M vs. M intcpts
# Z = beta(male vs female meta- analysis) - mean (of all male vs male meta-analysis estimates) / sd (of all male vs male meta-analysis estimates)
# Then calculate p-value from Z-score
# If p < 0.05 = male vs female intcpt is sig different to baseline male vs male intcpt


Z_MF_MM_intcpt <- (meta_intcpt_MF_across$b - MM_intcpt_means$mean_estimate) / MM_intcpt_means$mean_estimate_sd
p_MF_MM_intcpt <- 2 * pnorm(abs(Z_MF_MM_intcpt), mean = 0, sd = 1, lower.tail = F)



# Visualise
# Histogram of M vs M estimates and overlay M vs F intcpt
plot_distr_intcpt <- ggplot(meta_analyses_MM_intcpt_df, aes(x = estimate)) +
  geom_histogram(fill = "grey", colour = "black") +
  geom_vline(xintercept = meta_intcpt_MF_across$b[1,1], colour = "red") +
  geom_text(label = "Meta-analysis \nestimate of \nMale vs Female \nlinear regression \nintercept",
            x = 0.0003, y = 2500, colour = "red", hjust = 1) +
  scale_x_continuous("Meta-analysis estimates of Male vs Male linear regression intercept") +
  scale_y_continuous("Count") +
  theme_classic()

outfile = paste0(directory, "/male_vs_male_intcpts_distribution_overlaid_MF_intcpt.png")
ggsave(plot_distr_intcpt, file = outfile, width = 30, height = 30, unit = "cm")



# Create .txt output with results (incl slope and intercept)

# Define output file
outfile <- paste0(directory, "/meta_MF_vs_MM_linear_regression_summary.txt")
sink(outfile)

# SLOPE
cat("\n\nSlope from male vs female meta-analysis across cohorts vs distribution of slopes from male vs male meta-analyses:\n")

cat("\nDoes the observed M vs. F slope differ significantly from the distribution of M vs. M slopes?:\n")
cat("Z-score:", Z_MF_MM_slope, "\n")
cat("Two-tailed p-value:", p_MF_MM_slope, "\n")

if (p_MF_MM_slope < 0.05) {
  cat("Male vs female slope is significantly different to distribution of male vs male slopes.\n")
} else {
  cat("Male vs female slope is not significantly different to distribution of male vs male slopes.\n")
}


# INTERCEPT
cat("\n\nIntercept from male vs female meta-analysis across cohorts vs distribution of intercepts from male vs male meta-analyses:\n")

cat("\nDoes the observed M vs. F intercept differ significantly from the distribution of M vs. M intercepts?:\n")
cat("Z-score:", Z_MF_MM_intcpt, "\n")
cat("Two-tailed p-value:", p_MF_MM_intcpt, "\n")

if (p_MF_MM_intcpt < 0.05) {
  cat("Male vs female intercept is significantly different to distribution of male vs male intercepts.\n")
} else {
  cat("Male vs female intercept is not significantly different to distribution of male vs male intercepts.\n")
}


# Stop capturing output
sink()
## ----



## ---- CorrelationPlotRegressionFemalevsFemale
generate_plots_and_save_female_female_regression <- function(processed_dataframes) {
  # Define sample sizes for each cohort and sex
  sample_sizes <- list(
    AGDS = list(female = 17553),
    AllOfUs = list(female = 73454),
    Bionic = list(female = 37542),
    Blokland = list(female = 41937),
    GLAD = list(female = 20101),
    UKB = list(female = 99405)
  )
  
  cohort_names <- names(sample_sizes)
  
  # Generate all combinations of female cohorts (for female vs female comparison)
  combinations <- expand.grid(female_cohort1 = cohort_names, female_cohort2 = cohort_names)
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Create a data frame to store correlation results for female vs female
  regression_results <- data.frame(
    female_cohort1 = character(),
    female_cohort2 = character(),
    regression_slope = numeric(),
    regression_slope_se = numeric(),
    regression_slope_p = numeric(),
    regression_intercept = numeric(),
    regression_intercept_se = numeric(),
    regression_intercept_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each combination
  for (i in seq_len(nrow(combinations))) {
    female_cohort1 <- combinations$female_cohort1[i]
    female_cohort2 <- combinations$female_cohort2[i]
    
    # Extract relevant columns for the current combination
    if (female_cohort1 == female_cohort2) {
      data <- processed_dataframes %>%
        select(
          MARKER_build37,
          starts_with(paste0("BETA_", female_cohort1, "_female")),
          starts_with(paste0("SE_", female_cohort1, "_female")),
          starts_with(paste0("FREQA1_", female_cohort1, "_female"))
        ) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("FREQA1")) %>%
        mutate(
          BETA_female2 = BETA_female1,  # Copy BETA_female1 to BETA_female2
          SE_female2 = SE_female1,      # Copy SE_female1 to SE_female2
          FREQA1_female2 = FREQA1_female1  # Copy FREQA1_female1 to FREQA1_female2
        )
    } else {
      data <- processed_dataframes %>%
        select(
          MARKER_build37,
          starts_with(paste0("BETA_", female_cohort1, "_female")),
          starts_with(paste0("SE_", female_cohort1, "_female")),
          starts_with(paste0("BETA_", female_cohort2, "_female")),
          starts_with(paste0("SE_", female_cohort2, "_female")),
          starts_with(paste0("FREQA1_", female_cohort1, "_female")),
          starts_with(paste0("FREQA1_", female_cohort2, "_female"))
        ) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", female_cohort2, "_female"), "_female2", .), starts_with("BETA")) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", female_cohort2, "_female"), "_female2", .), starts_with("SE")) %>%
        rename_with(~ sub(paste0("_", female_cohort1, "_female"), "_female1", .), starts_with("FREQA1")) %>%
        rename_with(~ sub(paste0("_", female_cohort2, "_female"), "_female2", .), starts_with("FREQA1"))
    }
    
    
    # Standardize betas and SEs for both female cohorts
    female1_std <- calcu_std_b_se(
      z = data$BETA_female1 / data$SE_female1,
      p = data$FREQA1_female1,
      n = sample_sizes[[female_cohort1]]$female
    )
    female2_std <- calcu_std_b_se(
      z = data$BETA_female2 / data$SE_female2,
      p = data$FREQA1_female2,
      n = sample_sizes[[female_cohort2]]$female
    )
    
    # Create a new dataframe with standardized betas and standard errors
    standardized_data <- data.frame(
      MARKER_build37 = data$MARKER_build37,
      female1_std_b = female1_std$std_b_hat,
      female1_std_se = female1_std$std_se,
      female2_std_b = female2_std$std_b_hat,
      female2_std_se = female2_std$std_se
    )
    
    
    # Run linear regression of standardized female1 beta vs. standardized female2 beta
    lm_model <- lm(female2_std_b ~ female1_std_b, data = standardized_data)
    
    model_summary <- summary(lm_model)
    
    # Extract regression coefficients
    regression_slope <- model_summary$coefficients["female1_std_b", "Estimate"]
    regression_slope_se <- model_summary$coefficients["female1_std_b", "Std. Error"]
    regression_slope_p <- model_summary$coefficients["female1_std_b", "Pr(>|t|)"]
    regression_intercept <- model_summary$coefficients["(Intercept)", "Estimate"]
    regression_intercept_se <- model_summary$coefficients["(Intercept)", "Std. Error"]
    regression_intercept_p <- model_summary$coefficients["(Intercept)", "Pr(>|t|)"]
    
    # Add results to the regression_results data frame
    regression_results <- rbind(
      regression_results,
      data.frame(
        female_cohort1 = female_cohort1,
        female_cohort2 = female_cohort2,
        regression_slope = regression_slope,
        regression_slope_se = regression_slope_se,
        regression_slope_p = regression_slope_p,
        regression_intercept = regression_intercept,
        regression_intercept_se = regression_intercept_se,
        regression_intercept_p = regression_intercept_p
      )
    )
    
    # Calculate axis limits based on standardized values
    x_min <- min(standardized_data$female1_std_b - standardized_data$female1_std_se, na.rm = TRUE)
    x_max <- max(standardized_data$female1_std_b + standardized_data$female1_std_se, na.rm = TRUE)
    y_min <- min(standardized_data$female2_std_b - standardized_data$female2_std_se, na.rm = TRUE)
    y_max <- max(standardized_data$female2_std_b + standardized_data$female2_std_se, na.rm = TRUE)
    axis_limit <- max(abs(c(x_min, x_max, y_min, y_max)))
    label_x <- -0.9 * axis_limit
    label_y <- 0.9 * axis_limit
    
    # Generate the plot using standardized betas
    plot <- ggplot(standardized_data, aes(x = female1_std_b, y = female2_std_b)) +
      geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
      geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
      geom_linerange(aes(ymin = female2_std_b - female2_std_se, ymax = female2_std_b + female2_std_se), size = 0.4) +
      geom_linerangeh(aes(xmin = female1_std_b - female1_std_se, xmax = female1_std_b + female1_std_se), size = 0.4) +
      geom_point(alpha = 0.75) +
      geom_abline(colour = "grey40") +
      geom_smooth(method = "lm", se = FALSE, color = "grey80") +
      scale_x_continuous(paste("Standardized Female Beta (", female_cohort1, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      scale_y_continuous(paste("Standardized Female Beta (", female_cohort2, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      annotate("text", x = label_x, y = label_y,
               label = paste0(
                 "Intercept = ", round(regression_intercept, 3), " ± ", round(regression_intercept_se, 3),
                 ", P = ", format.pval(regression_intercept_p, digits = 3, scientific = TRUE), "\n",
                 "Beta = ", round(regression_slope, 3), " ± ", round(regression_slope_se, 3),
                 ", P = ", format.pval(regression_slope_p, digits = 3, scientific = TRUE)
               ),
               hjust = 0) +
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(female_cohort1, female_cohort2, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, regression = regression_results))
}


# Generate all plots
all_plots_regressions_female_v_female <- generate_plots_and_save_female_female_regression(processed_dataframes)

# Example: Display one plot
print(all_plots_regressions_female_v_female[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display correlations
all_plots_regressions_female_v_female$regression


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_regressions_female_v_female$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/female_vs_female_betas_linear_regression_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----


## ---- MetaAnalysisFemaleVsFemale

regression_female_v_female_df <- all_plots_regressions_female_v_female$regression %>%
  as_tibble() %>%
  mutate(cohort_label = paste(female_cohort1, "vs", female_cohort2))


# Create unique unordered pairs (15 independent pairs)
cohorts <- unique(c(regression_female_v_female_df$female_cohort1, regression_female_v_female_df$female_cohort2))
independent_pairs <- combn(cohorts, 2, simplify = FALSE)  # 15 unique unordered pairs

# Create both directions for each pair (A vs B and B vs A)
pairs_directions <- lapply(independent_pairs, function(pair) {
  list(paste(pair[1], pair[2], sep = " vs "), paste(pair[2], pair[1], sep = " vs "))
})

# Generate all possible sets of 15 pairs by selecting one direction per pair
# Each of the 15 pairs has 2 possible directions, so there are 2^15 = 32,768 sets
direction_combinations <- expand.grid(rep(list(c(1, 2)), 15))  # 32,768 rows, each a unique combination of directions (1st or 2nd option for each combo)

# Meta-analyse SLOPE
meta_results_slope_FF <- list()

for (i in seq_len(nrow(direction_combinations))) {
  selected_pairs <- sapply(seq_along(pairs_directions), function(j) {
    pairs_directions[[j]][[direction_combinations[i, j]]]  # Pick A->B or B->A for each pair
  })
  
  # Filter the regression dataset to include only selected pairs
  selected_data <- regression_female_v_female_df %>%
    filter(cohort_label %in% selected_pairs)
  
  # Run the meta-analysis
  meta_results_slope_FF[[paste("Meta_analysis_set", i, sep = "_")]] <- rma(
    yi = selected_data$regression_slope,
    sei = selected_data$regression_slope_se,
    weighted = TRUE,
    data = selected_data
  )
}

# Meta-analyse INTERCEPT
meta_results_intcpt_FF <- list()

for (i in seq_len(nrow(direction_combinations))) {
  selected_pairs <- sapply(seq_along(pairs_directions), function(j) {
    pairs_directions[[j]][[direction_combinations[i, j]]]  # Pick A->B or B->A for each pair
  })
  
  # Filter the regression dataset to include only selected pairs
  selected_data <- regression_female_v_female_df %>%
    filter(cohort_label %in% selected_pairs)
  
  # Run the meta-analysis
  meta_results_intcpt_FF[[paste("Meta_analysis_set", i, sep = "_")]] <- rma(
    yi = selected_data$regression_intercept,
    sei = selected_data$regression_intercept_se,
    weighted = TRUE,
    data = selected_data
  )
}

# View results
summary(meta_results_slope_FF[[1]])  # View first meta-analysis result SLOPE
summary(meta_results_intcpt_FF[[1]])  # View first meta-analysis result INTERCEPT

# Save results
# SLOPE
meta_analyses_FF_slope_df <- do.call(rbind, lapply(seq_along(meta_results_slope_FF), function(i) {
  meta <- meta_results_slope_FF[[i]]
  
  # Extract relevant values from the model summary
  data.frame(
    meta_set = paste("Meta_analysis_set", i, sep = "_"),
    tau2 = meta$tau2,   # Between-study variance
    tau = sqrt(meta$tau2),
    I2 = meta$I2,       # Percentage of variability due to heterogeneity
    H2 = meta$H2,       # Total variability / sampling variability
    Q = meta$QE,        # Cochran’s Q statistic
    Q_pval = meta$QEp,  # p-value for heterogeneity test
    estimate = meta$b,  # Meta-analysis estimated effect size
    se = meta$se,       # Standard error of the effect size
    sd = meta$se * sqrt(15),  # Standard deviation (SD = SE * sqrt(no. studies included = 15 in each))
    zval = meta$zval,   # Z-statistic
    pval = meta$pval,   # p-value for meta-analysis effect
    ci_lb = meta$ci.lb, # Lower bound of 95% CI
    ci_ub = meta$ci.ub  # Upper bound of 95% CI
  )
}))


outfile <- paste0(directory, "female_vs_female_meta_analyses_slope_all_32768_sets_of_15_linear_regressions_results.RData", sep = "")
save(meta_analyses_FF_slope_df, file = outfile)


# INTERCEPT
meta_analyses_FF_intcpt_df <- do.call(rbind, lapply(seq_along(meta_results_intcpt_FF), function(i) {
  meta <- meta_results_intcpt_FF[[i]]
  
  # Extract relevant values from the model summary
  data.frame(
    meta_set = paste("Meta_analysis_set", i, sep = "_"),
    tau2 = meta$tau2,   # Between-study variance
    tau = sqrt(meta$tau2),
    I2 = meta$I2,       # Percentage of variability due to heterogeneity
    H2 = meta$H2,       # Total variability / sampling variability
    Q = meta$QE,        # Cochran’s Q statistic
    Q_pval = meta$QEp,  # p-value for heterogeneity test
    estimate = meta$b,  # Meta-analysis estimated effect size
    se = meta$se,       # Standard error of the effect size
    sd = meta$se * sqrt(15),  # Standard deviation (SD = SE * sqrt(no. studies included = 15 in each))
    zval = meta$zval,   # Z-statistic
    pval = meta$pval,   # p-value for meta-analysis effect
    ci_lb = meta$ci.lb, # Lower bound of 95% CI
    ci_ub = meta$ci.ub  # Upper bound of 95% CI
  )
}))


outfile <- paste0(directory, "female_vs_female_meta_analyses_intcpt_all_32768_sets_of_15_linear_regressions_results.RData", sep = "")
save(meta_analyses_FF_intcpt_df, file = outfile)

# Calculate mean of all meta-analyses estimates, and sd of this mean estimate
FF_slope_means <- meta_analyses_FF_slope_df %>%
  summarise(mean_estimate = mean(estimate),
            mean_estimate_sd = sd(estimate),
            mean_estimate_se = mean_estimate_sd / sqrt(n()))

hist(meta_analyses_FF_slope_df$estimate)


FF_intcpt_means <- meta_analyses_FF_intcpt_df %>%
  summarise(mean_estimate = mean(estimate),
            mean_estimate_sd = sd(estimate),
            mean_estimate_se = mean_estimate_sd / sqrt(n()))

hist(meta_analyses_FF_intcpt_df$estimate)

# Define output file
outfile <- paste0(directory, "/meta_female_vs_female_across_cohorts_linear_regression_summary.txt")
sink(outfile)

## Slope
# Write the mean and sd of slope distributions
cat("\n\nResults from slope (beta value from regression of female beta vs female beta across cohorts) were meta-analysed in sets of 15 (32,768 meta-analyses):\n")

FF_slope_mean <- FF_slope_means %>% pull(mean_estimate)
FF_slope_sd <- FF_slope_means %>% pull(mean_estimate_sd)
FF_slope_se <- FF_slope_means %>% pull(mean_estimate_se)

cat("Mean slope from distribution of meta-anlyses", FF_slope_mean, "\n")
cat("SD of mean slope from distribution of meta-anlyses", FF_slope_sd, "\n")
cat("SE of mean slope from distribution of meta-anlyses", FF_slope_se, "\n")


## Intercept
cat("\n\nResults from intercept (beta value from regression of female beta vs female beta across cohorts) were meta-analysed in sets of 15 (32,768 meta-analyses):\n")

FF_intercept_mean <- FF_intcpt_means %>% pull(mean_estimate)
FF_intercept_sd <- FF_intcpt_means %>% pull(mean_estimate_sd)
FF_intercept_se <- FF_intcpt_means %>% pull(mean_estimate_se)

cat("Mean intercept from distribution of meta-anlyses", FF_intercept_mean, "\n")
cat("SD of mean intercept from distribution of meta-anlyses", FF_intercept_sd, "\n")
cat("SE of mean intercept from distribution of meta-anlyses", FF_intercept_se, "\n")

# Stop capturing output
sink()
## ----






## ---- CompareMFvsFF

## SLOPE

# Does the observed M vs. F slope differ significantly from the distribution of F vs. F slopes
# Z = beta(male vs female meta- analysis) - mean (of all female vs female meta-analysis estimates) / sd (of all female vs female meta-analysis estimates)
# Then calculate p-value from Z-score
# If p < 0.05 = male vs female slope is sig different to baseline female vs female slope


Z_MF_FF_slope <- (meta_slope_MF_across$b - FF_slope_means$mean_estimate) / FF_slope_means$mean_estimate_sd
p_MF_FF_slope <- 2 * pnorm(abs(Z_MF_FF_slope), mean = 0, sd = 1, lower.tail = F)



# Visualise
# Histogram of F vs F estimates and overlay M vs F slope
plot_distr_slope <- ggplot(meta_analyses_FF_slope_df, aes(x = estimate)) +
  geom_histogram(binwidth = 0.005, fill = "grey", colour = "black") +
  geom_vline(xintercept = meta_slope_MF_across$b[1,1], colour = "red") +
  geom_text(label = "Meta-analysis \nestimate of \nMale vs Female \nlinear regression \nslope",
            x = 0.3, y = 1500, colour = "red", hjust = 0, vjust = 1) +
  scale_x_continuous("Meta-analysis estimates of Female vs Female linear regression slope") +
  scale_y_continuous("Count") +
  theme_classic()

outfile = paste0(directory, "/female_vs_female_slopes_distribution_overlaid_MF_slope.png")
ggsave(plot_distr_slope, file = outfile, width = 30, height = 30, unit = "cm")




## INTERCEPT

# Does the observed M vs. F intcpt differ significantly from the distribution of F vs. F intcpts
# Z = beta(male vs female meta- analysis) - mean (of all female vs female meta-analysis estimates) / sd (of all female vs female meta-analysis estimates)
# Then calculate p-value from Z-score
# If p < 0.05 = male vs female intcpt is sig different to baseline female vs female intcpt


Z_MF_FF_intcpt <- (meta_intcpt_MF_across$b - FF_intcpt_means$mean_estimate) / FF_intcpt_means$mean_estimate_sd
p_MF_FF_intcpt <- 2 * pnorm(abs(Z_MF_FF_intcpt), mean = 0, sd = 1, lower.tail = F)



# Visualise
# Histogram of F vs F estimates and overlay M vs F intcpt
plot_distr_intcpt <- ggplot(meta_analyses_FF_intcpt_df, aes(x = estimate)) +
  geom_histogram(fill = "grey", colour = "black") +
  geom_vline(xintercept = meta_intcpt_MF_across$b[1,1], colour = "red") +
  geom_text(label = "Meta-analysis \nestimate of \nMale vs Female \nlinear regression \nintercept",
            x = 0.00035, y = 2500, colour = "red", hjust = 1) +
  scale_x_continuous("Meta-analysis estimates of Female vs Female linear regression intercept") +
  scale_y_continuous("Count") +
  theme_classic()

outfile = paste0(directory, "/female_vs_female_intcpts_distribution_overlaid_MF_intcpt.png")
ggsave(plot_distr_intcpt, file = outfile, width = 30, height = 30, unit = "cm")



# Create .txt output with results (incl slope and intercept)

# Define output file
outfile <- paste0(directory, "/meta_MF_vs_FF_linear_regression_summary.txt")
sink(outfile)

# SLOPE
cat("\n\nSlope from male vs female meta-analysis vs distribution of slopes from female vs female meta-analyses:\n")

cat("\nDoes the observed M vs. F slope differ significantly from the distribution of F vs. F slopes?:\n")
cat("Z-score:", Z_MF_FF_slope, "\n")
cat("Two-tailed p-value:", p_MF_FF_slope, "\n")

if (p_MF_FF_slope < 0.05) {
  cat("Male vs female slope is significantly different to distribution of female vs female slopes.\n")
} else {
  cat("Male vs female slope is not significantly different to distribution of female vs female slopes.\n")
}


# INTERCEPT
cat("\n\nIntercept from male vs female meta-analysis vs distribution of intercepts from female vs female meta-analyses:\n")

cat("\nDoes the observed M vs. F intercept differ significantly from the distribution of F vs. F intercepts?:\n")
cat("Z-score:", Z_MF_FF_intcpt, "\n")
cat("Two-tailed p-value:", p_MF_FF_intcpt, "\n")

if (p_MF_FF_intcpt < 0.05) {
  cat("Male vs female intercept is significantly different to distribution of female vs female intercepts.\n")
} else {
  cat("Male vs female intercept is not significantly different to distribution of female vs female intercepts.\n")
}


# Stop capturing output
sink()
## ----








## ---- MalevsFemaleGWASMASumstats
# Do male vs female correlation using M-A sumstats

# Restrict each sumstats to Adams lead independent genome-wide significant SNPs and make sure tested allele consistent for all sumstats
female_MA_dataset <- MA_female_df %>%
  mutate(A1 = toupper(Allele1),
         A2 = toupper(Allele2)) %>%
  inner_join(consistent_alleles, by = join_by("MarkerName" == "MARKER_build37")) %>%
  mutate(
    Flip = (A1 != Consistent_A1),
    BETA = ifelse(Flip, -Effect, Effect),
    A1 = ifelse(Flip, Consistent_A1, A1),
    A2 = ifelse(Flip, Consistent_A2, A2)
  ) %>%
  select(-Flip, -Consistent_A1, -Consistent_A2, -Effect) %>%
  rename_with(~paste0(., "_female"), -MarkerName)

male_MA_dataset <- MA_male_df %>%
  mutate(A1 = toupper(Allele1),
         A2 = toupper(Allele2)) %>%
  inner_join(consistent_alleles, by = join_by("MarkerName" == "MARKER_build37")) %>%
  mutate(
    Flip = (A1 != Consistent_A1),
    BETA = ifelse(Flip, -Effect, Effect),
    A1 = ifelse(Flip, Consistent_A1, A1),
    A2 = ifelse(Flip, Consistent_A2, A2)
  ) %>%
  select(-Flip, -Consistent_A1, -Consistent_A2, -Effect) %>%
  rename_with(~paste0(., "_male"), -MarkerName)


# Standardize betas and SEs for female and male cohorts
female_std <- calcu_std_b_se(
  z = female_MA_dataset$BETA_female / female_MA_dataset$StdErr_female,
  p = female_MA_dataset$Freq1_female,
  n = 289992
)
male_std <- calcu_std_b_se(
  z = male_MA_dataset$BETA_male / male_MA_dataset$StdErr_male,
  p = male_MA_dataset$Freq1_male,
  n = 196990
)

# Create a new dataframe with standardized betas and standard errors
standardized_data <- data.frame(
  MARKER_build37 = female_MA_dataset$MarkerName,
  female_std_b = female_std$std_b_hat,
  female_std_se = female_std$std_se,
  male_std_b = male_std$std_b_hat,
  male_std_se = male_std$std_se
)

save(standardized_data, file = paste0(directory, "Female_male_GWAS_meta_analysis_sumstats_standardised.RData"))


# Run linear regression of standardized male beta vs. standardized female beta

lm_MF_sumstats <- lm(male_std_b ~ female_std_b, data = standardized_data)

summary(lm_MF_sumstats)

MF_MA_df <- tidy(lm_MF_sumstats, conf.int = T, conf.level = 0.95)

slope <- MF_MA_df %>%
  filter(term == "female_std_b") %>%
  pull(estimate)

slope_se <- MF_MA_df %>%
  filter(term == "female_std_b") %>%
  pull(std.error)

intercept <- MF_MA_df %>%
  filter(term == "(Intercept)") %>%
  pull(estimate)

intercept_se <- MF_MA_df %>%
  filter(term == "(Intercept)") %>%
  pull(std.error)


# Define output file
outfile <- paste0(directory, "/male_vs_female_GWAS_MA_linear_regression_summary.txt")
sink(outfile)

## Slope
cat("\n\nSlope (from regression of male beta vs female beta using GWAS meta-analysis sumstats):\n")
cat("Slope:", slope, "\n")
cat("SE of lope:", slope_se, "\n")
print(summary(lm_MF_sumstats))

# estimate significantly different to 1?
hypothesized_value <- 1

Z <- (slope - hypothesized_value) / slope_se
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = TRUE)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 1):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 1 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 1 (p >= 0.05).\n")
}


## Intercept
cat("\n\nIntercept (from regression of male beta vs female beta using GWAS meta-analysis sumstats):\n")
cat("Intercept:", intercept, "\n")
cat("SE of lope:", intercept_se, "\n")
print(summary(lm_MF_sumstats))

# estimate significantly different to 0?
hypothesized_value <- 0

Z <- (intercept - hypothesized_value) / intercept_se
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = TRUE)

# Write Z-score and p-value
cat("\nHypothesis Testing (Null = estimate equals 0):\n")
cat("Z-score:", Z, "\n")
cat("Two-tailed p-value:", p, "\n")

if (p < 0.05) {
  cat("The estimate is significantly different to 1 (p < 0.05).\n")
} else {
  cat("The estimate is not significantly different to 1 (p >= 0.05).\n")
}


# Stop capturing output
sink()


# Calculate axis limits based on standardized values
x_min <- min(standardized_data$female_std_b - standardized_data$female_std_se, na.rm = TRUE)
x_max <- max(standardized_data$female_std_b + standardized_data$female_std_se, na.rm = TRUE)
y_min <- min(standardized_data$male_std_b - standardized_data$male_std_se, na.rm = TRUE)
y_max <- max(standardized_data$male_std_b + standardized_data$male_std_se, na.rm = TRUE)
axis_limit <- max(abs(c(x_min, x_max, y_min, y_max)))
label_x <- -0.9 * axis_limit
label_y <- 0.9 * axis_limit

# Generate the plot using standardized betas
plot <- ggplot(standardized_data, aes(x = female_std_b, y = male_std_b)) +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_linerange(aes(ymin = male_std_b - male_std_se, ymax = male_std_b + male_std_se), size = 0.4) +
  geom_linerangeh(aes(xmin = female_std_b - female_std_se, xmax = female_std_b + female_std_se), size = 0.4) +
  geom_point(alpha = 0.75) +
  geom_abline(colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "grey80") +
  scale_x_continuous(paste("Standardized Female Beta", sep = ""),
                     limits = c(-axis_limit, axis_limit)) +
  scale_y_continuous(paste("Standardized Male Beta", sep = ""),
                     limits = c(-axis_limit, axis_limit)) +
  annotate("text", x = label_x, y = label_y,
           label = paste("Slope =", round(slope, 3), "±", round(slope_se, 3),
                         "\nIntercept =", round(intercept, 3), "±", round(intercept_se, 3)),
           hjust = 0) +
  theme_classic()


plot

outfile = paste0(directory, "/male_vs_female_betas_linear_regression_MA_sumstats.png")
ggsave(plot, file = outfile, width = 20, height = 20, unit = "cm")
## ----



## ---- SummaryFigureSlope
# Male vs male data
MM_slope_means <- MM_slope_means %>%
  mutate(mean_estimate_ci_lower = mean_estimate - (1.96 * mean_estimate_se),
         mean_estimate_ci_upper = mean_estimate + (1.96 * mean_estimate_se),
         test = "MM",
         Z_1 = (mean_estimate - 1) / mean_estimate_se,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

# Female vs female data
FF_slope_means <- FF_slope_means %>%
  mutate(mean_estimate_ci_lower = mean_estimate - (1.96 * mean_estimate_se),
         mean_estimate_ci_upper = mean_estimate + (1.96 * mean_estimate_se),
         test = "FF",
         Z_1 = (mean_estimate - 1) / mean_estimate_se,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

# Male vs female data: across cohorts
MF_slope <- broom::tidy(meta_slope_MF_across,
                        conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_slope_across",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

MF_slope_across <- MF_slope %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high)


# Male vs female data: within same cohort
MF_slope <- broom::tidy(meta_slope_MF_within,
                        conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_slope_within",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

MF_slope_within <- MF_slope %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high)


# Male vs female data: Using GWAS sumstats
MF_MA <- MF_MA_df %>%
  filter(term == "female_std_b") %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high) %>%
  mutate(test = "MF_MA",
         Z_1 = (mean_estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

## Combine all data
all_data_slope_plot <- bind_rows(MM_slope_means, FF_slope_means, MF_slope_across, MF_slope_within, MF_MA)

all_data_slope_plot <- all_data_slope_plot %>%
  mutate(p_1_adj = p.adjust(p_1, method = "BH", n = 5))

save(all_data_slope_plot, file = paste0(directory, "all_data_for_plot_for_slopes.RData", sep = ""))
write.table(all_data_slope_plot, file = paste0(directory, "all_data_for_plot_for_slopes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


## Forest plot
plot_combined_slopes <- ggplot(all_data_slope_plot, aes(x = mean_estimate, y = test)) +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
  geom_vline(xintercept = 1,
             color = "grey80",
             linetype = "dashed") +
  scale_x_continuous("Linear regression slope",
                     limits = c(0, 1.1),
                     breaks = seq(0, 1, 0.2)) +
  scale_y_discrete("",
                   breaks = c("FF",
                              "MM",
                              "MF_slope_across",
                              "MF_slope_within",
                              "MF_MA"),
                   limits = c("FF",
                              "MM",
                              "MF_slope_across",
                              "MF_slope_within",
                              "MF_MA"),
                   labels = c("Female vs female (across cohorts)",
                              "Male vs male (across cohorts)",
                              "Male vs female (across cohorts)",
                              "Male vs female (within cohorts)",
                              "Male vs female (meta-analysed sumstats)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))


plot_combined_slopes

outfile <- paste(directory, "Forest_plot_all_comparisons_slopes.png", sep="")
ggsave(plot_combined_slopes, width = 20, height = 15, unit = "cm", file = outfile)

## ----




## ---- SummaryFigureIntercept

# Male vs male data
MM_intcpt_means <- MM_intcpt_means %>%
  mutate(mean_estimate_ci_lower = mean_estimate - (1.96 * mean_estimate_se),
         mean_estimate_ci_upper = mean_estimate + (1.96 * mean_estimate_se),
         test = "MM",
         Z_0 = (mean_estimate - 0) / mean_estimate_se,
         p_0 = 2 * pnorm(Z_0, mean = 0, sd = 1, lower.tail = F))

# Female vs female data
FF_intcpt_means <- FF_intcpt_means %>%
  mutate(mean_estimate_ci_lower = mean_estimate - (1.96 * mean_estimate_se),
         mean_estimate_ci_upper = mean_estimate + (1.96 * mean_estimate_se),
         test = "FF",
         Z_0 = (mean_estimate - 0) / mean_estimate_se,
         p_0 = 2 * pnorm(Z_0, mean = 0, sd = 1, lower.tail = F))

# Male vs female data: across cohorts
MF_intcpt <- broom::tidy(meta_intcpt_MF_across,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_intcpt_across",
         Z_0 = (estimate - 0) / std.error,
         p_0 = 2 * pnorm(Z_0, mean = 0, sd = 1, lower.tail = F))

MF_intcpt_across <- MF_intcpt %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high)


# Male vs female data: within same cohort
MF_intcpt <- broom::tidy(meta_intcpt_MF_within,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_intcpt_within",
         Z_0 = (estimate - 0) / std.error,
         p_0 = 2 * pnorm(Z_0, mean = 0, sd = 1, lower.tail = F))

MF_intcpt_within <- MF_intcpt %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high)


# Male vs female data: Using GWAS sumstats
MF_MA_intcpt <- MF_MA_df %>%
  filter(term == "(Intercept)") %>%
  rename(mean_estimate = estimate,
         mean_estimate_ci_lower = conf.low,
         mean_estimate_ci_upper = conf.high) %>%
  mutate(test = "MF_MA",
         Z_0 = (mean_estimate - 0) / std.error,
         p_0 = 2 * pnorm(Z_0, mean = 0, sd = 1, lower.tail = T))

## Combine all data
all_data_intcpt_plot <- bind_rows(MM_intcpt_means, FF_intcpt_means, MF_intcpt_across, MF_intcpt_within, MF_MA_intcpt)

all_data_intcpt_plot <- all_data_intcpt_plot %>%
  mutate(p_0_adj = p.adjust(p_0, method = "BH", n = 5))


save(all_data_intcpt_plot, file = paste0(directory, "all_data_for_plot_for_intcpts.RData", sep = ""))
write.table(all_data_intcpt_plot, file = paste0(directory, "all_data_for_plot_for_intcpts.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


## Forest plot
plot_combined_intcpts <- ggplot(all_data_intcpt_plot, aes(x = mean_estimate, y = test)) +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
  geom_vline(xintercept = 0,
             color = "grey80",
             linetype = "dashed") +
  scale_x_continuous("Linear regression intercept",
                     limits = c(-0.0004, 0.001),
                     breaks = seq(-0.0004, 0.001, 0.0002),
                     labels = function(x) ifelse(x == 0, "0", label_number(accuracy = 0.0001)(x))) +
  scale_y_discrete("",
                   breaks = c("FF",
                              "MM",
                              "MF_intcpt_across",
                              "MF_intcpt_within",
                              "MF_MA"),
                   limits = c("FF",
                              "MM",
                              "MF_intcpt_across",
                              "MF_intcpt_within",
                              "MF_MA"),
                   labels = c("Female vs female (across cohorts)",
                              "Male vs male (across cohorts)",
                              "Male vs female (across cohorts)",
                              "Male vs female (within cohorts)",
                              "Male vs female (meta-analysed sumstats)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))



plot_combined_intcpts

outfile <- paste(directory, "Forest_plot_all_comparisons_intcpts.png", sep="")
ggsave(plot_combined_intcpts, width = 20, height = 15, unit = "cm", file = outfile)

## ----


## ---- SummaryDataMFvsMMandMFvsFF
# Create summary data of comparing MF across cohorts to MM across cohorts and FF across cohorts, including p-value adjustment for multiple comparisons

comparison_slope_df <- tibble(
  comparison = c("MF_MM", "MF_FF"),
  measure = c("slope", "slope"),
  Z = c(Z_MF_MM_slope, Z_MF_FF_slope),
  p = c(p_MF_MM_slope, p_MF_FF_slope),
  p_adj = p.adjust(p, method = "BH", n = 2)
)

comparison_intcpt_df <- tibble(
  comparison = c("MF_MM", "MF_FF"),
  measure = c("intcpt", "intcpt"),
  Z = c(Z_MF_MM_intcpt, Z_MF_FF_intcpt),
  p = c(p_MF_MM_intcpt, p_MF_FF_intcpt),
  p_adj = p.adjust(p, method = "BH", n = 2)
)


write.table(comparison_slope_df, file = paste0(directory, "meta_MF_vs_MM_and_FF_linear_regression_slope_summary.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(comparison_intcpt_df, file = paste0(directory, "meta_MF_vs_MM_and_FF_linear_regression_intcpt_summary.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

## ----

