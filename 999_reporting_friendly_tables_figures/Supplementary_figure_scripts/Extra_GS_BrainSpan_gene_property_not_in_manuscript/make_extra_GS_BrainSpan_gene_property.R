# Create supplementary figure - Plots of developmental analysis run in FUMA using GxS sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Plot_developmental_GxS/"

Magma_dev_GxS = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/magma_exp_bs_dev_avg_log2RPKM.gsa.out"



### Packages ###
library(patchwork)
library(tidyverse)

### Load Data ###
Magma_dev_GxS_df <- read.table(Magma_dev_GxS, header = T, stringsAsFactors = F)



### Dataframe ###

# Combine the data from both GWAS' into a single dataframe
Magma_dev_GxS_df <- Magma_dev_GxS_df %>%
  mutate(VARIABLE = str_replace(VARIABLE, "_", " "))


# Reorder the VARIABLE factor based on the P values for the female dataset
female_order <- Magma_dev_GxS_df %>%
  arrange(P) %>%
  pull(VARIABLE)

# Apply the reordered factor levels to the whole dataset
Magma_dev_GxS_df <- Magma_dev_GxS_df %>%
  mutate(VARIABLE = factor(VARIABLE, levels = female_order))


### Plot ###

plot <- ggplot(data = Magma_dev_GxS_df, aes(x = VARIABLE, y = -log10(P), fill = Dataset)) +
  geom_col(position = position_dodge(),
           fill = "#20A387FF") +
  geom_hline(yintercept = -log10(0.05/11),
             colour = "grey30",
             linetype = 'dotted') +
  labs(x = "",
       y = expression(-log[10](p-value))) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        legend.position = "top",
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))

plot



outfile <- paste(directory, "Bar_gene_property_developmental_analysis_GxS.png", sep="")
ggsave(plot, width = 24, height = 10, unit = "cm", file = outfile)
