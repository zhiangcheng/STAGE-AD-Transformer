# Create supplementary figure - Manhattan plot of gene-based tests run in FUMA using GxS fulldc sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Manhattan_gene_based_GxS_fulldc/"

Magma_genes_GxS = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/magma.genes.out"


### Packages ###
library(ggrepel)      # To repel labels on plot
library(tidyverse)

### Load Data ###
Magma_genes_GxS_df <- read.table(Magma_genes_GxS, header = T, stringsAsFactors = F)


# Create new column to use for labelling significant genes
# And change chr X to 23
Magma_genes_GxS_df <- Magma_genes_GxS_df %>%
  mutate(gene = case_when(
    P < 2.53e-06 ~ SYMBOL
  ),
  CHR = if_else(CHR == 'X', '23', CHR),
  CHR = as.numeric(CHR))

# Create cumulative BP position for x axis
cum_BP <- Magma_genes_GxS_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(START)) %>%
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

Man_plot_df <- Magma_genes_GxS_df %>%
  inner_join(cum_BP, by = "CHR") %>%
  mutate(bp_cum = START + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

# Set max y axis value
ylim <- Man_plot_df %>%
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)


# Set significance threshold
sig <- 2.53e-06





Man_plot <- ggplot(Man_plot_df,
                     aes(x = bp_cum,
                         y = -log10(P),
                         color = as_factor(CHR)))+
  geom_hline(aes(yintercept = -log10(sig)),
             color = "grey30",
             linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(-ylim, ylim, 1),
                     labels = abs(seq(-ylim, ylim, 1))) +
  scale_color_manual(values = rep(c("#20A387FF", "#95D840FF"),
                                  unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome",
       y = expression(-log[10](p-value))) +
  geom_text_repel(data = subset(Man_plot_df, !is.na(gene)),
                  aes(label = gene),
                  size = 3, family = "Calibri", color = "#191d1f",
                  box.padding = 0.25, point.padding = 0, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "#191d1f", segment.size = 0.3) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "#191d1f", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "#191d1f", margin = margin(0,10,0,0)),
        axis.line = element_line(colour = "#191d1f"),
        axis.ticks = element_line(colour = "#191d1f"))


Man_plot

outfile <- paste(directory, "Manhattan_plot_gene_based_GxS_fulldc.png", sep="")
ggsave(Man_plot, width = 19, height = 10, unit = "cm", file = outfile)



