# This script tests if female and male betas are different
# Does correlation of female-female, male-male and female-male betas for SNPs known to be associated with depression in each cohort
# Then meta-analyses each of these R across all the cohorts

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
library(MBESS)      # SE of R2
library(psych)      # for correlation
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


# Function for calculating pearsons correlation
calc_pearson_cor <- function(x,y){
  result <- corr.test(x, y, use = "complete", method = "pearson")
  return(result$r)
}

calc_pearson_cor_se <- function(x,y){
  result <- corr.test(x, y, use = "complete", method = "pearson")
  return(result$se)
}

## ----


## ---- CorrelationPlotRegressionMalevsFemale
generate_plots_and_save_cor <- function(processed_dataframes) {
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
  cor_results <- data.frame(
    female_cohort = character(),
    male_cohort = character(),
    pearson_r = numeric(),
    pearson_r_se = numeric(),
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
    
    # Calculate pearsons correlation between standardized betas
    pearson_r <- calc_pearson_cor(
      x = standardized_data$female_std_b,
      y = standardized_data$male_std_b
    )
    
    pearson_r_se <- calc_pearson_cor_se(
      x = standardized_data$female_std_b,
      y = standardized_data$male_std_b
    )
    
    
    # Add results to the cor_results data frame
    cor_results <- rbind(
      cor_results,
      data.frame(
        female_cohort = female_cohort,
        male_cohort = male_cohort,
        pearson_r = pearson_r,
        pearson_r_se = pearson_r_se
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
               label = bquote(R == .(round(pearson_r, 3)) ~ "±" ~ .(round(pearson_r_se, 3))),
               hjust = 0) +
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(female_cohort, male_cohort, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, correlations = cor_results))
}


# Generate all plots
all_plots_cor_male_v_female <- generate_plots_and_save_cor(processed_dataframes)

# Example: Display one plot
print(all_plots_cor_male_v_female[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display regression results
all_plots_cor_male_v_female$correlations


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_cor_male_v_female$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/male_vs_female_betas_correlation_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----


## ---- MetaAnalysisMaleVsFemale
# Meta-analyse (inverse variance weighted) the r for the M-F regression in each of the 6 cohorts
# Uses random effects (heterogeneity between cohorts)

cor_male_v_female_df <- all_plots_cor_male_v_female$correlations %>%
  as_tibble() %>%
  mutate(cohort_label = paste(female_cohort, "vs", male_cohort))

save(cor_male_v_female_df, file = paste(directory, "Correlation_data_male_vs_female.RData", sep = ""))

## Across cohorts (all comparisons apart from diagonal = 30)

cor_male_v_female_across_df <- cor_male_v_female_df %>%
  filter(female_cohort != male_cohort)

# Meta-analysis of male beta vs female beta across cohorts r2 value
meta_r_MF_across <- rma(yi = pearson_r,
                        sei = pearson_r_se,
                        weighted = T,
                        data = cor_male_v_female_across_df)


summary(meta_r_MF_across)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_r_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_r_MF_across, slab = cor_male_v_female_across_df$cohort_label, xlab = "R of male beta vs female beta")
dev.off()


# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_summary.txt")
sink(outfile)

# Write the summary
cat("\n\nSummary of meta-analysis of R (from correlation of male beta vs female beta across cohorts):\n")
print(summary(meta_r_MF_across))

# estimate significantly different to 1?
estimate <- summary(meta_r_MF_across)[["beta"]]
SE <- summary(meta_r_MF_across)[["se"]]
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

# Stop capturing output
sink()




## Within cohorts (all diagonal comparisons = 6)

cor_male_v_female_within_df <- cor_male_v_female_df %>%
  filter(female_cohort == male_cohort)

# Meta-analysis of male beta vs female R2 within same cohort
meta_r_MF_within <- rma(yi = pearson_r,
                        sei = pearson_r_se,
                        weighted = T,
                        data = cor_male_v_female_within_df)


summary(meta_r_MF_within)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_r_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_r_MF_within, slab = cor_male_v_female_within_df$cohort_label, xlab = "R of male beta vs female beta")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_summary.txt")
sink(outfile)

# Write the summary of
cat("\n\nSummary of meta-analysis of R (from correlation of male beta vs female beta within the same cohort):\n")
print(summary(meta_r_MF_within))

# estimate significantly different to 1?
estimate <- summary(meta_r_MF_within)[["beta"]]
SE <- summary(meta_r_MF_within)[["se"]]
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


# Stop capturing output
sink()
## ----






## ---- CorrelationPlotRegressionMalevsMale
generate_plots_and_save_male_male_correlation <- function(processed_dataframes) {
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
  correlation_results <- data.frame(
    male_cohort1 = character(),
    male_cohort2 = character(),
    pearson_r = numeric(),
    pearson_r_se = numeric(),
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
    
    
    # Calculate pearsons correlation between standardized betas
    pearson_r <- calc_pearson_cor(
      x = standardized_data$male1_std_b,
      y = standardized_data$male2_std_b
    )
    
    pearson_r_se <- calc_pearson_cor_se(
      x = standardized_data$male1_std_b,
      y = standardized_data$male2_std_b
    )
    
    
    # Add results to the correlation_results data frame
    correlation_results <- rbind(
      correlation_results,
      data.frame(
        male_cohort1 = male_cohort1,
        male_cohort2 = male_cohort2,
        pearson_r = pearson_r,
        pearson_r_se = pearson_r_se
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
               label = bquote(R == .(round(pearson_r, 3)) ~ "±" ~ .(round(pearson_r_se, 3))),
               hjust = 0) +
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(male_cohort1, male_cohort2, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, correlation = correlation_results))
}


# Generate all plots
all_plots_correlations_male_v_male <- generate_plots_and_save_male_male_correlation(processed_dataframes)

# Example: Display one plot
print(all_plots_correlations_male_v_male[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display correlations
all_plots_correlations_male_v_male$correlation


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_correlations_male_v_male$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/male_vs_male_betas_correlation_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----


## ---- MetaAnalysisMaleVsMale
# Meta-analyse (inverse variance weighted) the correlations for the M-M correlation in each of the 6 cohorts
# Uses random effects (heterogeneity between cohorts)

## Across cohorts, above diagonal (15 comparisons - as x and y axes swapped give same R)
correlation_male_v_male_df <- all_plots_correlations_male_v_male$correlation %>%
  as_tibble() %>%
  mutate(cohort_label = paste(male_cohort1, "vs", male_cohort2))

save(correlation_male_v_male_df, file = paste(directory, "Correlation_data_male_vs_male.RData", sep = ""))

# Top half of 6x6 plot grid is when the x axis cohort (male_cohort1) comes alphabetically before y axis cohort (male_cohort2)
correlation_male_v_male_across_df <- correlation_male_v_male_df %>%
  filter(as.character(male_cohort1) < as.character(male_cohort2))


# Meta-analysis of male beta vs male beta correlation
meta_r_MM_across <- rma(yi = pearson_r,
                        sei = pearson_r_se,
                        weighted = T,
                        data = correlation_male_v_male_across_df)


summary(meta_r_MM_across)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_male_across_r_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_r_MM_across, slab = correlation_male_v_male_across_df$cohort_label, xlab = "R of male beta vs male beta")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_male_vs_male_across_cohort_summary.txt")
sink(outfile)


# Write the summary
cat("\n\nSummary of meta-analysis of R (from correlation of male beta vs male beta across cohorts):\n")
print(summary(meta_r_MM_across))

# estimate significantly different to 1?
estimate <- summary(meta_r_MM_across)[["beta"]]
SE <- summary(meta_r_MM_across)[["se"]]
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


# Stop capturing output
sink()
## ----







## ---- CorrelationPlotRegressionFemalevsFemale
generate_plots_and_save_female_female_correlation <- function(processed_dataframes) {
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
  correlation_results <- data.frame(
    female_cohort1 = character(),
    female_cohort2 = character(),
    pearson_r = numeric(),
    pearson_r_se = numeric(),
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
    
    # Calculate pearsons correlation between standardized betas
    pearson_r <- calc_pearson_cor(
      x = standardized_data$female1_std_b,
      y = standardized_data$female2_std_b
    )
    
    pearson_r_se <- calc_pearson_cor_se(
      x = standardized_data$female1_std_b,
      y = standardized_data$female2_std_b
    )
    
    
    # Add results to the correlation_results data frame
    correlation_results <- rbind(
      correlation_results,
      data.frame(
        female_cohort1 = female_cohort1,
        female_cohort2 = female_cohort2,
        pearson_r = pearson_r,
        pearson_r_se = pearson_r_se
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
      scale_x_continuous(paste("Standardized Male Beta (", female_cohort1, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      scale_y_continuous(paste("Standardized Male Beta (", female_cohort2, ")", sep = ""),
                         limits = c(-axis_limit, axis_limit)) +
      annotate("text", x = label_x, y = label_y,
               label = bquote(R == .(round(pearson_r, 3)) ~ "±" ~ .(round(pearson_r_se, 3))),
               hjust = 0) +
      
      theme_classic()
    
    # Save plot in the list
    plot_list[[paste(female_cohort1, female_cohort2, sep = "_vs_")]] <- plot
  }
  
  # Return both plots and correlation results
  return(list(plots = plot_list, correlation = correlation_results))
}


# Generate all plots
all_plots_correlations_female_v_female <- generate_plots_and_save_female_female_correlation(processed_dataframes)

# Example: Display one plot
print(all_plots_correlations_female_v_female[["plots"]][["AGDS_vs_Bionic"]])

# Example: Display correlations
all_plots_correlations_female_v_female$correlation


# Generate the 6x6 grid of plots
plot_grid <- marrangeGrob(
  grobs = all_plots_correlations_female_v_female$plots, # List of plots
  nrow = 6, # Number of rows in the grid
  ncol = 6 # Number of columns in the grid
)


outfile = paste0(directory, "/female_vs_female_betas_correlation_36_plots.png")
ggsave(plot_grid, file = outfile, width = 80, height = 80, unit = "cm")
## ----



## ---- MetaAnalysisFemaleVsFemale
# Meta-analyse (inverse variance weighted) the correlations for the F-F in each of the 6 cohorts
# Uses random effects (heterogeneity between cohorts)

## Across cohorts, above diagonal (15 comparisons, as flipped x and y axes give same R)
correlation_female_v_female_df <- all_plots_correlations_female_v_female$correlation %>%
  as_tibble() %>%
  mutate(cohort_label = paste(female_cohort1, "vs", female_cohort2))

save(correlation_female_v_female_df, file = paste(directory, "Correlation_data_female_vs_female.RData", sep = ""))

# Top half of 6x6 plot grid is when the x axis cohort (female_cohort1) comes alphabetically before y axis cohort (female_cohort2)
correlation_female_v_female_across_df <- correlation_female_v_female_df %>%
  filter(as.character(female_cohort1) < as.character(female_cohort2))


# Meta-analysis of female beta vs female beta
meta_r_FF_across <- rma(yi = pearson_r,
                        sei = pearson_r_se,
                        weighted = T,
                        data = correlation_female_v_female_across_df)


summary(meta_r_FF_across)

# Forest plot
outfile <- paste0(directory, "/meta_female_vs_female_across_r_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_r_FF_across, slab = correlation_female_v_female_across_df$cohort_label, xlab = "R of female beta vs female beta")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_female_vs_female_across_cohort_summary.txt")
sink(outfile)


# Write the summary
cat("\n\nSummary of meta-analysis of R (from correlation of female beta vs female beta across cohorts):\n")
print(summary(meta_r_FF_across))

# estimate significantly different to 1?
estimate <- summary(meta_r_FF_across)[["beta"]]
SE <- summary(meta_r_FF_across)[["se"]]
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


# Run correlation of standardized male beta vs. standardized female beta
pearson_r <- calc_pearson_cor(
  x = standardized_data$female_std_b,
  y = standardized_data$male_std_b
)

pearson_r_se <- calc_pearson_cor_se(
  x = standardized_data$female_std_b,
  y = standardized_data$male_std_b
)


result <- corr.test(standardized_data$female_std_b,
                    standardized_data$male_std_b,
                    use = "complete",
                    method = "pearson")

r <- result$r
se <- result$se
lower_ci <- as.numeric(result$ci$lower)
upper_ci <- as.numeric(result$ci$upper)

MF_MA_df <- tibble(
  estimate = as.numeric(r),
  std.error = as.numeric(se),
  conf.low = lower_ci,
  conf.high = upper_ci,
  test = "MF_MA"
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
  scale_x_continuous(paste("Standardized Female Beta", sep = ""),
                     limits = c(-axis_limit, axis_limit)) +
  scale_y_continuous(paste("Standardized Male Beta", sep = ""),
                     limits = c(-axis_limit, axis_limit)) +
  annotate("text", x = label_x, y = label_y,
           label = bquote(R == .(round(pearson_r, 3)) ~ "±" ~ .(round(pearson_r_se, 3))),
           hjust = 0) +
  theme_classic()


plot

outfile = paste0(directory, "/male_vs_female_betas_correlation_MA_sumstats.png")
ggsave(plot, file = outfile, width = 20, height = 20, unit = "cm")
## ----




## ---- SummaryFigure

# Male vs male data
MM_r <- broom::tidy(meta_r_MM_across,
                    conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MM_r_across",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

MM_r_across_I2 <- confint(meta_r_MM_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MM_r_across <- bind_cols(MM_r, MM_r_across_I2)


# Female vs female data
FF_r <- broom::tidy(meta_r_FF_across,
                    conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "FF_r_across",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

FF_r_across_I2 <- confint(meta_r_FF_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

FF_r_across <- bind_cols(FF_r, FF_r_across_I2)


# Male vs female data: across cohorts
MF_r <- broom::tidy(meta_r_MF_across,
                    conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_r_across",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

MF_r_across_I2 <- confint(meta_r_MF_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MF_r_across <- bind_cols(MF_r, MF_r_across_I2)


# Male vs female data: within same cohort
MF_r <- broom::tidy(meta_r_MF_within,
                    conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_r_within",
         Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))

MF_r_within_I2 <- confint(meta_r_MF_within, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MF_r_within <- bind_cols(MF_r, MF_r_within_I2)


# Male vs female data: Using GWAS M-A susmtats

MF_MA_df <- MF_MA_df %>%
  mutate(Z_1 = (estimate - 1) / std.error,
         p_1 = 2 * pnorm(Z_1, mean = 0, sd = 1, lower.tail = T))


## Combine all data
all_data <- bind_rows(MM_r_across, FF_r_across, MF_r_across, MF_r_within, MF_MA_df)
all_data <- all_data %>%
  mutate(p_1_adj = p.adjust(p_1, method = "BH", n = 5))

save(all_data, file = paste0(directory, "all_data_for_plot.RData", sep = ""))
write.table(all_data, file = paste0(directory, "all_data_for_plot.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

## Forest plot
plot_combined <- ggplot(all_data, aes(x = estimate, y = test)) +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1,
             color = "grey80",
             linetype = "dashed") +
  geom_text(aes(x = 1.1, y = test,
                label = ifelse(test == "MF_MA", NA, paste0(round(`estimate_I^2(%)`, 0), "%")))) +
  scale_x_continuous("Correlation of female and male betas",
                     limits = c(0, 1.1),
                     breaks = seq(0, 1, 0.2)) +
  scale_y_discrete("",
                   breaks = c("FF_r_across",
                              "MM_r_across",
                              "MF_r_across",
                              "MF_r_within",
                              "MF_MA"),
                   limits = c("FF_r_across",
                              "MM_r_across",
                              "MF_r_across",
                              "MF_r_within",
                              "MF_MA"),
                   labels = c("Female vs female (across cohorts)",
                              "Male vs male (across cohorts)",
                              "Male vs female (across cohorts)",
                              "Male vs female (within cohorts)",
                              "Male vs female (meta-analysed sumstats)")) +
  annotate("text", x = 1.1, y = 5.4,
           label = bquote(I^2),
           hjust = 0) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))


plot_combined

outfile <- paste(directory, "Forest_plot_all_comparisons.png", sep="")
ggsave(plot_combined, width = 20, height = 15, unit = "cm", file = outfile)

## ----
