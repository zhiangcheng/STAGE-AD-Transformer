# Create supplementary figure with QQ plot, lambda and LDSC results from GxS GWAS

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(patchwork)    # To create panelled figure
library(tidyverse)

### Load Data ###
GxS_fulldc <- "/path/03_Metal/GxS/fulldc/Metaanalysis_AD_GxS_fulldc_AllCohorts_QCed_rsID.txt"
GxS_nodc <- "/path/03_Metal/GxS/nodc/Metaanalysis_AD_GxS_nodc_AllCohorts_QCed_rsID.txt"

GxS_fulldc_df <- read.table(GxS_fulldc, header = TRUE, stringsAsFactors = F)
GxS_nodc_df <- read.table(GxS_nodc, header = TRUE, stringsAsFactors = F)

### QQ Plot fulldc ###
# calculate lambda
alpha_fulldc <- median(qchisq(1-GxS_fulldc_df$P,1))/qchisq(0.5,1)

# specify confidence interval
ci <- 0.95
# no. SNPs
n_snps_fulldc <- nrow(GxS_fulldc_df)

# Create dataframe with all the data for the QQ Plot
QQ_plot_df_fulldc <- tibble(
  observed = -log10(sort(GxS_fulldc_df$P)),
  expected = -log10(ppoints(n_snps_fulldc)),
  clower = -log10(qbeta(
    p = (1 - ci) / 2,
    shape1 = seq(n_snps_fulldc),
    shape2 = rev(seq(n_snps_fulldc))
  )),
  cupper = -log10(qbeta(
    p = (1 + ci) / 2,
    shape1 = seq(n_snps_fulldc),
    shape2 = rev(seq(n_snps_fulldc))
  ))
)

#Plot
QQ_plot_fulldc <- ggplot(QQ_plot_df_fulldc, aes(x = expected, y = observed)) +
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
       caption = paste("lambda =", signif(alpha_fulldc, digits = 4),
                       "\n LDSC intercept = 1.0162 (95% CI: 1.0052 - 1.0272)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.caption = element_text(size = 12))



### QQ Plot nodc ###
# calculate lambda
alpha_nodc <- median(qchisq(1-GWAS_nodc_df$P,1))/qchisq(0.5,1)

# specify confidence interval
ci <- 0.95
# no. SNPs
n_snps_nodc <- nrow(GWAS_nodc_df)

# Create dataframe with all the data for the QQ Plot
QQ_plot_df_nodc <- tibble(
  observed = -log10(sort(GWAS_nodc_df$P)),
  expected = -log10(ppoints(n_snps_nodc)),
  clower = -log10(qbeta(
    p = (1 - ci) / 2,
    shape1 = seq(n_snps_nodc),
    shape2 = rev(seq(n_snps_nodc))
  )),
  cupper = -log10(qbeta(
    p = (1 + ci) / 2,
    shape1 = seq(n_snps_nodc),
    shape2 = rev(seq(n_snps_nodc))
  ))
)

#Plot
QQ_plot_nodc <- ggplot(QQ_plot_df_nodc, aes(x = expected, y = observed)) +
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
       caption = paste("lambda =", signif(alpha_nodc, digits = 4), 
                       "\n LDSC intercept = 1.0162 (95% CI: 1.0052 - 1.0272)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.caption = element_text(size = 12))


## Panelled Figure ##

figure <- QQ_plot_fulldc + QQ_plot_nodc +
  plot_annotation(tag_levels = 'A')

outfile <- paste("/path/10_Publication_figures/Supplementary_figures/QQ_plots_lambda_LDSCintcpt/QQ_plots_GxS.png", sep="")
ggsave(figure, width = 20, height = 10, unit = "cm", file = outfile)


