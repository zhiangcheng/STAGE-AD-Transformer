# Create supplementary fugure - Manhattan plots of gene-based tests run in FUMA using sex-stratified GWAS sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Miami_gene_based_female_vs_male/"

Magma_genes_female = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/magma.genes.out"

Magma_genes_male = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Males_genomewidesig_positionalmap10kb_MHCexcl/magma.genes.out"

### Packages ###
library(data.table)   # To load data
library(ggh4x)        # To create different scales on ggplot facets
library(ggrepel)      # To repel lables on plot
library(tidyverse)

### Load Data ###
Magma_genes_female_df <- read.table(Magma_genes_female, header = T, stringsAsFactors = F)

Magma_genes_male_df <- read.table(Magma_genes_male, header = T, stringsAsFactors = F)

# Combine the data from both GWAS' into a single dataframe
combined_df <- bind_rows(
  mutate(Magma_genes_female_df, Dataset = "female"),
  mutate(Magma_genes_male_df, Dataset = "male")
)

# Create new column to use for labelling significant genes
# And change chr X to 23
combined_df <- combined_df %>% 
  mutate(gene = case_when(
    P < 2.53e-06 ~ SYMBOL
  ),
  CHR = if_else(CHR == 'X', '23', CHR),
  CHR = as.numeric(CHR))

# Create cumulative BP position for x axis
cum_BP <- combined_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(START)) %>%
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

Man_plot_df <- combined_df %>%
  inner_join(cum_BP, by = "CHR") %>%
  mutate(bp_cum = START + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

# Set max y axis value
ylim <- Man_plot_df %>%
  filter(Dataset == "female") %>%
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "male") %>%
  filter(P == min(P)) %>%
  mutate(ylim = (abs(floor(log10(P))) + 1)) %>%
  pull(ylim)

# Use the max ylim so the two plots are comparable
chosen_ylim <- max(ylim, ylim2)

# Set significance threshold
sig <- 2.53e-06


# Plot
facet_labels <- c(
  female = paste("Females"),
  male = paste("Males")
)


Miami_plot <- ggplot(Man_plot_df,
                     aes(x = bp_cum,
                         y = ifelse(Dataset == "female",
                                    -log10(P), log10(P)),
                         color = as_factor(CHR)))+
  geom_hline(aes(yintercept = ifelse(Dataset == "female",
                                     -log10(sig), log10(sig))),
             color = "grey30",
             linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(-chosen_ylim, chosen_ylim, 1),
                     labels = abs(seq(-chosen_ylim, chosen_ylim, 1))) +
  scale_color_manual(values = rep(c("#20A387FF", "#95D840FF"),
                                  unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome",
       y = expression(-log[10](p-value)),
       color = "Dataset") +
  geom_text_repel(data = subset(Man_plot_df, !is.na(gene)),
                  aes(label = gene),
                  size = 3, family = "Calibri", color = "#191d1f",
                  box.padding = 0.25, point.padding = 0, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "#191d1f", segment.size = 0.3) +
  facet_wrap(~ Dataset, scale = "free_y", ncol = 1, strip.position = "right",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "#191d1f", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "#191d1f", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "#191d1f"),
        axis.line = element_line(colour = "#191d1f"),
        axis.ticks = element_line(colour = "#191d1f"))

# To set y axis scales for the two facets with the same chosen_ylim so have same y axis making it easier to compare
position_scales <- list(
  scale_y_continuous(limits = c(0, chosen_ylim),
                     breaks = seq(0, chosen_ylim, 2),
                     labels = abs(seq(0, chosen_ylim, 2))),
  scale_y_continuous(limits = c(-chosen_ylim, 0),
                     breaks = seq(-chosen_ylim, 0, 2),
                     labels = abs(seq(-chosen_ylim, 0, 2))))

Miami_plot <- Miami_plot + facetted_pos_scales(y = position_scales)

Miami_plot

outfile <- paste(directory, "Miami_plot_gene_based_sex_stratified.png", sep="")
ggsave(Miami_plot, width = 19, height = 10, unit = "cm", file = outfile)






