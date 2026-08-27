# Create supplementary figure with QQ plot, lambda and LDSC results from sex-stratified GWAS

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(patchwork)    # To create panelled figure
library(tidyverse)

### Load Data ###
GWAS_female <- "/path/03_Metal/Females/Metaanalysis_AD_female_AllCohorts_QCed_rsID.txt"
GWAS_male <- "/path/03_Metal/Males/Metaanalysis_AD_male_AllCohorts_QCed_rsID.txt"

GWAS_female_df <- read.table(GWAS_female, header = TRUE, stringsAsFactors = F)
GWAS_male_df <- read.table(GWAS_male, header = TRUE, stringsAsFactors = F)

### QQ Plot Female ###
# calculate lambda
alpha_female <- median(qchisq(1-GWAS_female_df$P,1))/qchisq(0.5,1)

# specify confidence interval
ci <- 0.95
# no. SNPs
n_snps_female <- nrow(GWAS_female_df)

# Create dataframe with all the data for the QQ Plot
QQ_plot_df_female <- tibble(
  observed = -log10(sort(GWAS_female_df$P)),
  expected = -log10(ppoints(n_snps_female)),
  clower = -log10(qbeta(
    p = (1 - ci) / 2,
    shape1 = seq(n_snps_female),
    shape2 = rev(seq(n_snps_female))
  )),
  cupper = -log10(qbeta(
    p = (1 + ci) / 2,
    shape1 = seq(n_snps_female),
    shape2 = rev(seq(n_snps_female))
  ))
)

#Plot
QQ_plot_female <- ggplot(QQ_plot_df_female, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower),
              fill = "grey30", alpha = 0.5) +
  geom_point() +
  geom_segment(data = . %>% filter(expected == max(expected)),
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5,
               color = "grey30", lineend = "round") +
  scale_x_continuous(lim = c(0, 8.5),
                     breaks = seq(0, 8, 1)) +
  scale_y_continuous(lim = c(0, 10), 
                     breaks = seq(0, 10, 1)) +
  labs(x = expression(Expected  -log[10](p-value)),
       y = expression(Observed  -log[10](p-value)),
       caption = paste("lambda =", signif(alpha_female, digits = 4),
                       "\n LDSC intercept = 1.0064 (95% CI: 0.9898 - 1.0230)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.caption = element_text(size = 12))



### QQ Plot Male ###
# calculate lambda
alpha_male <- median(qchisq(1-GWAS_male_df$P,1))/qchisq(0.5,1)

# specify confidence interval
ci <- 0.95
# no. SNPs
n_snps_male <- nrow(GWAS_male_df)

# Create dataframe with all the data for the QQ Plot
QQ_plot_df_male <- tibble(
  observed = -log10(sort(GWAS_male_df$P)),
  expected = -log10(ppoints(n_snps_male)),
  clower = -log10(qbeta(
    p = (1 - ci) / 2,
    shape1 = seq(n_snps_male),
    shape2 = rev(seq(n_snps_male))
  )),
  cupper = -log10(qbeta(
    p = (1 + ci) / 2,
    shape1 = seq(n_snps_male),
    shape2 = rev(seq(n_snps_male))
  ))
)

#Plot
QQ_plot_male <- ggplot(QQ_plot_df_male, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower),
              fill = "grey30", alpha = 0.5) +
  geom_point() +
  geom_segment(data = . %>% filter(expected == max(expected)),
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5,
               color = "grey30", lineend = "round") +
  scale_x_continuous(lim = c(0, 8.5),
                     breaks = seq(0, 8, 1)) +
  scale_y_continuous(lim = c(0, 10), 
                     breaks = seq(0, 10, 1)) +
  labs(x = expression(Expected  -log[10](p-value)),
       y = expression(Observed  -log[10](p-value)),
       caption = paste("lambda =", signif(alpha_male, digits = 4), 
                       "\n LDSC intercept = 1.0014 (95% CI: 0.9864 - 1.0164)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.caption = element_text(size = 12))


## Panelled Figure ##

figure <- QQ_plot_female + QQ_plot_male +
  plot_annotation(tag_levels = 'A')

outfile <- paste("/path/10_Publication_figures/Supplementary_figures/QQ_plots_lambda_LDSCintcpt/QQ_plots_female_male.png", sep="")
ggsave(figure, width = 20, height = 10, unit = "cm", file = outfile)


