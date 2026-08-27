# Ran LDSC rg for all female and male cohort pairwise comparisons
# Now meta-analyse M-F, F-F and M-M

directory = "/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/"

rg_results = "/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/ldsc_rg_MDD_all_female_male_cohorts_results.txt"

rg_MA_results <- "/path/04_LDSC/SNPrg/Females/rg_femaleMDD_correlation_table.txt"

### Packages ###

library(metafor)      # meta-analysis
library(broom)   # create tibble from rma output
library(tidyverse)



### Load Data ###

rg_results_df <- read.table(rg_results, header = T)

rg_results_df <- rg_results_df %>%
  mutate(
    cohort1 = sub("_.*", "", p1),
    sex1 = sub(".*_", "", p1),
    cohort2 = sub("_.*", "", p2),
    sex2 = sub(".*_", "", p2)
  ) %>%
  mutate(cohort_label = paste(cohort1, sex1, "vs", cohort2, sex2))

rg_MA_results_df <- read.table(rg_MA_results, header = T)

rg_MA_results_df <- rg_MA_results_df %>%
  filter(p1 == "Allcohorts_female_MDD" & p2 == "Allcohorts_male_MDD")




### Male vs Female: Across Cohorts ###

# All male vs female comparisons that are across cohorts (no male vs female comparisons within the same cohort)
# male vs female and female vs male when the same cohorts are involved will have the same rg so only include one of these

rg_results_male_female_across <- rg_results_df %>%
  filter(sex1 == "female" & sex2 == "male",
         cohort1 != cohort2)

# Meta-analysis of rg
meta_rg_MF_across <- rma(yi = rg,
                         sei = se,
                         weighted = T,
                         data = rg_results_male_female_across)


summary(meta_rg_MF_across)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_rg_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_rg_MF_across, slab = rg_results_male_female_across$cohort_label, xlab = "Rg of male vs female across cohorts")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_across_cohorts_summary.txt")
sink(outfile)

cat("\n\nSummary of meta-analysis of rg male vs female across cohorts:\n")
print(summary(meta_rg_MF_across))

# estimate significantly different to 1?
estimate <- summary(meta_rg_MF_across)[["beta"]]
SE <- summary(meta_rg_MF_across)[["se"]]
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



### Male vs Female: Within Same Cohort ###

# All male vs female comparisons that are within the same cohort
# male vs female and female vs male when the same cohorts are involved will have the same rg so only include one of these

rg_results_male_female_within <- rg_results_df %>%
  filter(sex1 == "female" & sex2 == "male",
         cohort1 == cohort2)

# Meta-analysis of rg
meta_rg_MF_within <- rma(yi = rg,
                         sei = se,
                         weighted = T,
                         data = rg_results_male_female_within)


summary(meta_rg_MF_within)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_rg_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_rg_MF_within, slab = rg_results_male_female_within$cohort_label, xlab = "Rg of male vs female within same cohort")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_male_vs_female_within_cohort_summary.txt")
sink(outfile)

cat("\n\nSummary of meta-analysis of rg male vs female within same cohort:\n")
print(summary(meta_rg_MF_within))

# estimate significantly different to 1?
estimate <- summary(meta_rg_MF_within)[["beta"]]
SE <- summary(meta_rg_MF_within)[["se"]]
hypothesized_value <- 1

Z <- (estimate - hypothesized_value) / SE
p <- 2 * pnorm(Z, mean = 0, sd = 1, lower.tail = FALSE)

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



### Female vs female: Across Cohorts ###

# All female vs female comparisons that are across cohorts (no female vs female comparisons within the same cohort as would just be 1 as comparing the same sumstats with itself)
# female cohort 1 vs female cohort 2 and female cohort 2 vs female cohort 1 will have the same rg so only include one of these

rg_results_female_female_across <- rg_results_df %>%
  filter(sex1 == "female" & sex2 == "female",
         cohort1 != cohort2)

# Meta-analysis of rg
meta_rg_FF_across <- rma(yi = rg,
                         sei = se,
                         weighted = T,
                         data = rg_results_female_female_across)


summary(meta_rg_FF_across)

# Forest plot
outfile <- paste0(directory, "/meta_female_vs_female_across_cohorts_rg_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_rg_FF_across, slab = rg_results_female_female_across$cohort_label, xlab = "Rg of female vs female across cohorts")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_female_vs_female_across_cohorts_summary.txt")
sink(outfile)

cat("\n\nSummary of meta-analysis of rg female vs female across cohorts:\n")
print(summary(meta_rg_FF_across))

# estimate significantly different to 1?
estimate <- summary(meta_rg_FF_across)[["beta"]]
SE <- summary(meta_rg_FF_across)[["se"]]
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


### Male vs Male: Across Cohorts ###

# All male vs male comparisons that are across cohorts (no male vs male comparisons within the same cohort as would just be 1 as comparing the same sumstats with itself)
# male cohort 1 vs male cohort 2 and male cohort 2 vs male cohort 1 will have the same rg so only include one of these

rg_results_male_male_across <- rg_results_df %>%
  filter(sex1 == "male" & sex2 == "male",
         cohort1 != cohort2)

# Meta-analysis of rg
meta_rg_MM_across <- rma(yi = rg,
                         sei = se,
                         weighted = T,
                         data = rg_results_male_male_across)


summary(meta_rg_MM_across)

# Forest plot
outfile <- paste0(directory, "/meta_male_vs_male_across_cohorts_rg_forest_plot.png", sep = "")
png(outfile, width = 30, height = 25, units = "cm", res = 300)
forest(meta_rg_MM_across, slab = rg_results_male_male_across$cohort_label, xlab = "Rg of male vs male across cohorts")
dev.off()

# Define output file
outfile <- paste0(directory, "/meta_male_vs_male_across_cohorts_summary.txt")
sink(outfile)

cat("\n\nSummary of meta-analysis of rg male vs male across cohorts:\n")
print(summary(meta_rg_MM_across))

# estimate significantly different to 1?
estimate <- summary(meta_rg_MM_across)[["beta"]]
SE <- summary(meta_rg_MM_across)[["se"]]
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




### Summary Figure ###

# Male vs male data
MM_across <- broom::tidy(meta_rg_MM_across,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MM_across")

MM_across_I2 <- confint(meta_rg_MM_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MM_across <- bind_cols(MM_across, MM_across_I2)


# female vs female data
FF_across <- broom::tidy(meta_rg_FF_across,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "FF_across")

FF_across_I2 <- confint(meta_rg_FF_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

FF_across <- bind_cols(FF_across, FF_across_I2)


# Male vs female data
MF_across <- broom::tidy(meta_rg_MF_across,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_across")

MF_across_I2 <- confint(meta_rg_MF_across, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MF_across <- bind_cols(MF_across, MF_across_I2)



MF_within <- broom::tidy(meta_rg_MF_within,
                         conf.int = T, conf.level = 0.95,) %>%
  mutate(test = "MF_within")

MF_within_I2 <- confint(meta_rg_MF_within, fixed = F, random = T, level = 95) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = rowname, values_from = c(estimate, ci.lb, ci.ub))

MF_within <- bind_cols(MF_within, MF_within_I2)


# Male vs female using GWAS meta-analysis sumstats data
MF_MA <- rg_MA_results_df %>%
  rename(estimate = rg,
         std.error = se,
         statistic = z,
         p.value = p) %>%
  mutate(test = "MF_MA",
         conf.low = estimate - (1.96 * std.error),
         conf.high = estimate + (1.96 * std.error)) %>%
  select(estimate, std.error, statistic, p.value, conf.low, conf.high, test)

all_data <- bind_rows(MM_across, FF_across, MF_across, MF_within, MF_MA)

all_data <- all_data %>%
  mutate(p.value_adj = p.adjust(p.value, method = "BH", n = 5))

save(all_data, file = paste0(directory, "all_data_for_plot.RData", sep = ""))


plot_estimate <- ggplot(all_data, aes(x = estimate, y = test)) +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1,
             color = "grey80",
             linetype = "dashed") +
  geom_text(aes(x = 1.2, y = test,
                label = ifelse(test == "MF_MA", NA, paste0(round(`estimate_I^2(%)`, 0), "%")))) +
  scale_x_continuous("Genetic Correlation",
                     limits = c(0, 1.25),
                     breaks = seq(0, 1.2, 0.2)) +
  scale_y_discrete("",
                   breaks = rev(c("MF_MA",
                                  "MF_within",
                                  "MF_across",
                                  "MM_across",
                                  "FF_across")),
                   limits = rev(c("MF_MA",
                                  "MF_within",
                                  "MF_across",
                                  "MM_across",
                                  "FF_across")),
                   labels = rev(c("Male vs Female (meta-analysed sumstats)",
                                  "Male vs Female (within cohorts)",
                                  "Male vs Female (across cohorts)",
                                  "Male vs Male (across cohorts)",
                                  "Female vs Female (across cohorts)"))) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))


plot_estimate

outfile <- paste(directory, "Metaanalysis_rg_MF_MM_FF_across_cohorts.png", sep="")
ggsave(plot_estimate, width = 20, height = 15, unit = "cm", file = outfile)



