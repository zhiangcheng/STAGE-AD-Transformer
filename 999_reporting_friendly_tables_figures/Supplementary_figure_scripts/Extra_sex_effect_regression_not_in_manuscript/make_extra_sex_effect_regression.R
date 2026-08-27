# Create publication figure of forest plots for male beta vs female beta linear regression slope and intercept

directory = "/path/10_Publication_figures/Supplementary_figures/Betas_female_vs_male_linear_regression_forest_plots/"



### Packages ###
library(patchwork)
library(scales)      # change scientific notation on plots
library(tidyverse)


### Load Data ###
# Linear regression data
load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot_for_slopes.RData")

load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot_for_intcpts.RData")

# Load all cohorts linear regression data - so can include distribution in graphs

# Slope F-F
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/female_vs_female_meta_analyses_slope_all_32768_sets_of_15_linear_regressions_results.RData")
# Intcpt F-F
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/female_vs_female_meta_analyses_intcpt_all_32768_sets_of_15_linear_regressions_results.RData")

# Slope M-M
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/male_vs_male_meta_analyses_slope_all_32768_sets_of_15_linear_regressions_results.RData")
# Intcpt M-M
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/male_vs_male_meta_analyses_intcpt_all_32768_sets_of_15_linear_regressions_results.RData")

# Slope and Intcpt F-M 
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/male_vs_female_linear_regression_results_all_cohort_combos.RData")



# Slope df for plot
slope_FF <- meta_analyses_FF_slope_df %>% 
  mutate(test = "FF")

slope_MM <- meta_analyses_MM_slope_df %>% 
  mutate(test = "MM")

slope_within_sex <- bind_rows(
  slope_FF,
  slope_MM
)

slope_MF <- regression_male_v_female_df %>% 
  mutate(test = case_when(
    female_cohort != male_cohort ~ "MF_slope_across",
    female_cohort == male_cohort ~ "MF_slope_within"
  ))

# Intcpt df for plot
intcpt_FF <- meta_analyses_FF_intcpt_df %>% 
  mutate(test = "FF")

intcpt_MM <- meta_analyses_MM_intcpt_df %>% 
  mutate(test = "MM")

intcpt_within_sex <- bind_rows(
  intcpt_FF,
  intcpt_MM
)

intcpt_MF <- regression_male_v_female_df %>% 
  mutate(test = case_when(
    female_cohort != male_cohort ~ "MF_intcpt_across",
    female_cohort == male_cohort ~ "MF_intcpt_within"
  ))


## Meta-Analysis of Female vs Male Betas Plot: Linear regression slope ###

Slope_plot <- ggplot(all_data_slope_plot, aes(x = mean_estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_violin(data = slope_within_sex, aes(x = estimate, y = test, colour = test, fill = test),
              inherit.aes = F,
              show.legend = F) +
  geom_jitter(data = slope_MF, aes(x = regression_slope, y = test, colour = test),
              alpha = 0.2, height = 0.1,
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
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
                   labels = c("Female vs female \n(across cohorts)",
                              "Male vs male \n(across cohorts)",
                              "Male vs female \n(across cohorts)",
                              "Male vs female \n(within cohorts)",
                              "Male vs female \n(meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      breaks = c("FF",
                                 "MM",
                                 "MF_slope_across",
                                 "MF_slope_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  scale_x_continuous("Linear regression slope") +
  annotate("text", x = 0.9, y = "MF_MA", label = "*", size = 6) +
  annotate("text", x = 0.59, y = "MF_slope_within", label = "*", size = 6) +
  annotate("text", x = 0.52, y = "MF_slope_across", label = "*", size = 6) +
  annotate("text", x = 0.5, y = "MM", label = "*", size = 6) +
  annotate("text", x = 0.6, y = "FF", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

Slope_plot



### Meta-Analysis of Female vs Male Betas Plot: Linear regression intercept ###
dynamic_labels <- function(x) {
  sapply(x, function(val) {
    if (val == 0) {
      "0"
    } else {
      # Count number of decimal places
      dp <- nchar(sub("^[^.]*\\.?([^0]*?)0*$", "\\1", format(val, scientific = FALSE)))
      acc <- 10^(-dp)
      label_number(accuracy = acc)(val)
    }
  })
}

Intcpt_plot <- ggplot(all_data_intcpt_plot, aes(x = mean_estimate, y = test, colour = test)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_violin(data = intcpt_within_sex, aes(x = estimate, y = test, colour = test, fill = test),
              inherit.aes = F,
              show.legend = F) +
  geom_jitter(data = intcpt_MF, aes(x = regression_intercept, y = test, colour = test),
              alpha = 0.2, height = 0.1,
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
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
  scale_colour_manual("Comparison",
                      breaks = c("FF",
                                 "MM",
                                 "MF_intcpt_across",
                                 "MF_intcpt_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  scale_x_continuous("Linear regression intercept",
                     limits = c(-0.0016, 0.003),
                     breaks = seq(-0.002, 0.003, 0.001),
                     labels = function(x) ifelse(x == 0, "0", label_number(accuracy = 0.001)(x))) +
  annotate("text", x = 0.001, y = "MF_intcpt_across", label = "*", size = 6) +
  annotate("text", x = 0.0015, y = "MM", label = "*", size = 6) +
  annotate("text", x = 0.0009, y = "FF", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

Intcpt_plot


### Create Panelled Figure ###

figure <- Slope_plot | Intcpt_plot+
  plot_annotation(tag_levels = 'A')


figure

outfile <- paste(directory, "Male_vs_female_betas_linear_regression.png", sep="")
ggsave(figure, width = 21, height = 12, unit = "cm", file = outfile)
