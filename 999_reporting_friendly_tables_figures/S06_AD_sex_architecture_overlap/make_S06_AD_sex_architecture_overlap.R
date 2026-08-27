# Create publication figure of tests comparing males and females
# SBayesS (heritability, polygenicity, selection paramter), rg, correlation, polygenic overlap (rg, mixer)

directory = "/path/10_Publication_figures/02_Compare_sexes/"

# All MCMC samples so can plot distribution

SBayes_directory = "/path/12_SBayesS_h2/"

MA_females_mcmc <- paste0(SBayes_directory, "Metaanalysis_AD_female_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

MA_males_mcmc <- paste0(SBayes_directory, "Metaanalysis_AD_male_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

# All cohort rg results so can plot distribution

rg_results_cohorts = "/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/ldsc_rg_AD_all_female_male_cohorts_results.txt"


### Packages ###
library(patchwork)
library(ggforce)      # ggplot circles
library(showtext)
showtext_auto()
library(tidyverse)


### Load Data ###

# heritability data
load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_MA <- h2_liability_varying_K


# Polygenicity data
load("/path/12_SBayesS_h2/Pi_data_MA_cohorts.RData")


# Selection parameter data
load("/path/12_SBayesS_h2/S_data_MA_cohorts.RData")


# Genetic correlation data
load("/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/all_data_for_plot.RData")
all_data_rg <- all_data

# Pearson correlation data
load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot.RData")
all_data_R <- all_data

# Load all mcmc reports

MA_females_mcmc_df <- read.table(MA_females_mcmc, header = TRUE, stringsAsFactors = FALSE)
MA_males_mcmc_df <- read.table(MA_males_mcmc, header = TRUE, stringsAsFactors = FALSE)

MA_females_mcmc_df <- MA_females_mcmc_df %>%
  mutate(sex = "female")

MA_males_mcmc_df <- MA_males_mcmc_df %>%
  mutate(sex = "male")

MA_all_mcmc <- bind_rows(MA_females_mcmc_df, MA_males_mcmc_df)

# Remove burn-in of 5,000 (as thinning of 10 used = 500 rows)
MA_all_mcmc <- MA_all_mcmc %>%
  slice(501:n())

# Load all cohorts rg data
rg_results_cohorts_df <- read.table(rg_results_cohorts, header = T)

rg_results_cohorts_df <- rg_results_cohorts_df %>%
  mutate(
    cohort1 = sub("_.*", "", p1),
    sex1 = sub(".*_", "", p1),
    cohort2 = sub("_.*", "", p2),
    sex2 = sub(".*_", "", p2)
  ) %>%
  mutate(cohort_label = paste(cohort1, sex1, "vs", cohort2, sex2))

# Load all cohorts correlation data
load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/Correlation_data_female_vs_female.RData")

load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/Correlation_data_male_vs_male.RData")

load(file = "/path/09_Effect_sizes_plots/Adams_SNPs/Correlation_data_male_vs_female.RData")

# Make into one df
correlation_male_v_male_df <- correlation_male_v_male_df %>%
  mutate(test = case_when(
    as.character(male_cohort1) < as.character(male_cohort2) ~ "MM_r_across"
  )) %>%
  filter(test == "MM_r_across")

correlation_female_v_female_df <- correlation_female_v_female_df %>%
  mutate(test = case_when(
    as.character(female_cohort1) < as.character(female_cohort2) ~ "FF_r_across"
  )) %>%
  filter(test == "FF_r_across")

cor_male_v_female_df <- cor_male_v_female_df %>%
  mutate(test = case_when(
    female_cohort != male_cohort ~ "MF_r_across",
    female_cohort == male_cohort ~ "MF_r_within"
  ))

cor_results_cohorts_df <- bind_rows(correlation_male_v_male_df, correlation_female_v_female_df, cor_male_v_female_df)

### Calculate h2 on liability scale ###

# Function to convert h2 to liability scale
# Lee et al 2011 AJHG 's method to convert the heritability estimate and standard error at the observed scale to those at the liability scale obs is the estimate (or SE) at the observed scale K is the population prevalence P is the sample prevalence

mapToLiabilityScale = function(obs, K, P){
  z = dnorm(qnorm(1-K))
  lia = obs * (K*(1-K)/z^2) * (K*(1-K)/(P*(1-P)))
  return(lia)
}

# Convert all hsq values in MCMC samples to liability scale
MA_all_mcmc <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.20, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.10, P = 0.3290))
  )

### SNP-based heritability Plot: Liability scale ###

# Using results from SBayesS on GWAS sex stratified meta-analysis results

h2_liability <- h2_liability_varying_K_MA %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(scale = "Liability") %>%
  select(-K)


h2_plot_sexes <- ggplot(h2_liability, aes(x = sex, y = mean_hsq * 100, colour = sex)) +
  geom_violin(data = MA_all_mcmc, aes(x = sex, y = hsq_l * 100, colour = sex, fill = sex),
              inherit.aes = F,
              show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100),
                 show.legend = F) +
  geom_point(data = data.frame(sex = 'across', mean_hsq = NA)) +
  geom_linerange(data = data.frame(sex = 'across', dummy_x = 'male', ymin = -999, ymax = -999),
                 aes(x = dummy_x, ymin = ymin, ymax = ymax, colour = sex),
                 inherit.aes = FALSE, show.legend = F) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(8, 12.5)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_manual(
    values = c('female' = '#FDE725FF',
               'male' = '#440154FF',
               'across' = '#20A387FF'),
    breaks = c("male", "female", "across"),
    labels = c('male' = "Male",
               'female' = "Female",
               'across' = "Across sex")
  ) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  annotate("text", x = 1.5, y = 12.5, label = "P(F > M) = 100%", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"))

h2_plot_sexes


### Polygenicity Plot ###

pi_plot <- ggplot(all_data_pi %>% filter(Type == "MA"), aes(x = sex, y = mean_pi, colour = sex)) +
  geom_violin(data = MA_all_mcmc, aes(x = sex, y = Pi, colour = sex, fill = sex),
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_y_continuous(paste("Polygenicity (\u03c0)"),
                     limits = c(0,0.031)) +
  annotate("text", x = 1.5, y = 0.03, label = "P(F > M) = 100%", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"),
        legend.position = "none")

pi_plot


### Selection Parameter Plot ###

S_plot <- ggplot(all_data_S %>% filter(Type == "MA"), aes(x = sex, y = mean_S, colour = sex)) +
  geom_violin(data = MA_all_mcmc, aes(x = sex, y = S, colour = sex, fill = sex),
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_y_continuous(expression(paste("Selection parameter (" * italic(S) * ")")),
                     limits = c(-0.35, 0.23),
                     breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),
                     labels = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2)) +
  annotate("text", x = 1.5, y = 0.23, label = "P(F > M) = 78%", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"),
        legend.position = "none")

S_plot



### Female - male genetic correlation plot ###

# For all cohort data add 'test' column so aligns in plot
rg_results_cohorts_df <- rg_results_cohorts_df %>%
  mutate(test = case_when(
    sex1 == "female" & sex2 == "male" & cohort1 != cohort2 ~ "MF_across",
    sex1 == "female" & sex2 == "male" & cohort1 == cohort2 ~ "MF_within",
    sex1 == "female" & sex2 == "female" & cohort1 != cohort2 ~ "FF_across",
    sex1 == "male" & sex2 == "male" & cohort1 != cohort2 ~ "MM_across",
  ))

rg_plot <- ggplot(all_data_rg, aes(x = estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_jitter(data = rg_results_cohorts_df, aes(x = rg, y = test, colour = test),
              alpha = 0.2, height = 0.1,
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_y_discrete("",
                   breaks = c("FF_across",
                              "MM_across",
                              "MF_across",
                              "MF_within",
                              "MF_MA"),
                   limits = c("FF_across",
                              "MM_across",
                              "MF_across",
                              "MF_within",
                              "MF_MA"),
                   labels = c("Female vs female \n(across cohorts)",
                              "Male vs male \n(across cohorts)",
                              "Male vs female \n(across cohorts)",
                              "Male vs female \n(within cohorts)",
                              "Male vs female \n(meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      values = c("FF_across" = "#FDE725FF",
                                 "MM_across" = "#440154FF",
                                 "MF_across" = "#20A387FF",
                                 "MF_within" = "#20A387FF",
                                 "MF_MA" = "#20A387FF"),
                      breaks = c("FF_across", "MM_across", "MF_across"),
                      labels = c("Female", "Male", "Across sex")) +
  guides(color = "none") +
  scale_x_continuous(expression(italic(r[g])),
                     limits = c(0, 1.55),
                     breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5),
                     labels = c("0", "0.25", "0.5", "0.75", "1", "1.25", "1.5")) +
  annotate("text", x = 0.97, y = "MF_MA",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.82, y = "MF_across",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.87, y = "FF_across",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"),
        legend.position = "none")

rg_plot


### Meta-Analysis of Female vs Male Betas Plot: Correlation (R) ###

R_plot <- ggplot(all_data_R, aes(x = estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_jitter(data = cor_results_cohorts_df, aes(x = pearson_r, y = test, colour = test),
              alpha = 0.2, height = 0.1,
              inherit.aes = F, show.legend = F) +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
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
  scale_colour_manual("Comparison",
                      breaks = c("FF_r_across",
                                 "MM_r_across",
                                 "MF_r_across",
                                 "MF_r_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  guides(color = "none") +
  scale_x_continuous("Correlation",
                     limits = c(0, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  annotate("text", x = 0.87, y = "MF_MA", vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.55, y = "MF_r_within",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.4, y = "MF_r_across",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.37, y = "MM_r_across",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  annotate("text", x = 0.5, y = "FF_r_across",  vjust = -0.6, label = "*", size = 7 / 2.845276) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"),
        legend.position = "none")

R_plot


### Mixer Venn ###

# Data: sizes of each group
female <- 6133
male <- 0
overlap <- 7111

sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Calculate circle radii proportional to the square root of sizes
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate offset to shift Male circle to the right so its right edge touches Female's right edge
offset <- r_female - r_male

# Create circle data with adjusted position for Male circle
circle_data <- data.frame(
  x = c(1, 1 + offset),  # Offset the Male circle to the right
  y = c(1, 1),           # Keep y-coordinates the same for female and male circle
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram
mixer_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  guides(fill = "none") +
  coord_fixed() +
  # labs(title = "Number of causal variants") +
  annotate("text", x = circle_data$x[1] - 50, y = circle_data$y[1],
           label = female, size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text", x = circle_data$x[2], y = circle_data$y[2],
           label = overlap, size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn




### gwas-pw Venn ###

# Data: sizes of each group
female <- 3
male <- 0
overlap <- 42

sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Calculate circle radii proportional to the square root of sizes
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate offset to shift Male circle to the right so its right edge touches Female's right edge
offset <- r_female - r_male

# Create circle data with adjusted position for Male circle
circle_data <- data.frame(
  x = c(1, 1 + offset),  # Offset the Male circle to the right
  y = c(1, 1),           # Keep y-coordinates the same for female and male circle
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram
gwaspw_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  guides(fill = "none") +
  coord_fixed() +
  # labs(title = "Number of genomic regions") +
  annotate("text", x = circle_data$x[1] - 4.5, y = circle_data$y[1],
           label = female, size = 7 / 2.845276, color = "black", family = "Arial", fontface = "bold") +
  annotate("text", x = circle_data$x[2], y = circle_data$y[2],
           label = overlap, size = 7 / 2.845276, color = "white", family = "Arial", fontface = "bold") +
  geom_segment(aes(x = circle_data$x[1] - 4.2, xend = circle_data$x[1] - 3.65,
                   y = circle_data$y[1], yend = circle_data$y[1]),
               color = "black", size = 0.5) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

gwaspw_venn


### Create Panelled Figure ###

figure <- (guide_area() /
             ((h2_plot_sexes | pi_plot | S_plot ) + plot_layout(widths = c(1, 1, 1))) /
             ((rg_plot | R_plot) + plot_layout(widths = c(1, 1))) /
             ((mixer_venn | gwaspw_venn) + plot_layout(widths = c(1, 1)))) +
  plot_layout(heights = c(0.2, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = "top")

figure

outfile <- paste(directory, "Compare_sexes_h2_pi_S_rg_R_mixer_gwaspw.eps", sep="")
ggsave(figure, width = 18, height = 18.5, units = "cm", file = outfile, device = cairo_ps, dpi = 300)

