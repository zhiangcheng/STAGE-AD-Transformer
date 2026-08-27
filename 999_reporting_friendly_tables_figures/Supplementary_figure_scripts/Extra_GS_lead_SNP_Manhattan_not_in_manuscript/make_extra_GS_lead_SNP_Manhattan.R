directory = "/path/10_Publication_figures/Supplementary_figures/Manhattan_GxS_Adams_SNPs/"

GxS_GWAS_Adams_sig <- "/path/08_Sensitivity_analyses/Adams2024_topSNPs/Metaanalysis_AD_GxS_fulldc_AllCohorts_QCed_rsID_Adams_sig_SNPs.txt"

#########################################################################

### Packages ###
library(tidyverse)

### Load Data ###
GxS_GWAS_Adams_sig_df <- read.table(GxS_GWAS_Adams_sig, header = T, stringsAsFactors = F)


### Manhattan Plot ###

# Create cumulative BP position for x axis
cum_BP <- GxS_GWAS_Adams_sig_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

Man_plot_df <- GxS_GWAS_Adams_sig_df %>%
  inner_join(cum_BP, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

# Set significance threshold
sig <- 0.05/697


## Man Plot: P-value GxE ##

# Set max y axis value
ylimP <- 5



# Plot
Man_plot <- ggplot(Man_plot_df,
                   aes(x = bp_cum, y = -log10(P), color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig),
             color = "grey30",
             linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, ylimP),
                     breaks = seq(0, ylimP, 1),
                     labels = abs(seq(0, ylimP, 1))) +
  scale_color_manual(values = rep(c("#20A387FF", "#95D840FF"),
                                  unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome",
       y = expression(-log[10](p-value))) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))


Man_plot

outfile <- paste(directory, "Manhattan_plot_GxS_fulldc_Adams_sig_SNPs.png", sep="")
ggsave(Man_plot, width = 19, height = 10, unit = "cm", file = outfile)

