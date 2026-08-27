# Create supplementary fugure - Plots of developmental analysis run in FUMA using sex-stratified GWAS sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Plot_developmental_female_male/"

Magma_dev_female = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_bs_dev_avg_log2RPKM.gsa.out"

Magma_dev_male = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Males_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_bs_dev_avg_log2RPKM.gsa.out"


### Packages ###
library(patchwork)
library(tidyverse)

### Load Data ###
Magma_dev_female_df <- read.table(Magma_dev_female, header = T, stringsAsFactors = F)

Magma_dev_male_df <- read.table(Magma_dev_male, header = T, stringsAsFactors = F)


### Dataframe ###

# Combine the data from both GWAS' into a single dataframe
combined_df <- bind_rows(
  mutate(Magma_dev_female_df, Dataset = "female"),
  mutate(Magma_dev_male_df, Dataset = "male")
) %>%
  mutate(VARIABLE = str_replace(VARIABLE, "_", " "))


# Reorder the VARIABLE factor based on the P values for the female dataset
female_order <- combined_df %>%
  filter(Dataset == "female") %>%
  arrange(P) %>%
  pull(VARIABLE)

# Apply the reordered factor levels to the whole dataset
combined_df <- combined_df %>%
  mutate(VARIABLE = factor(VARIABLE, levels = female_order))


### Plot ###

plot <- ggplot(data = combined_df, aes(x = VARIABLE, y = -log10(P), fill = Dataset)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d("Sex", begin = 1, end = 0,
                       labels = c("Females", "Males")) +
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



outfile <- paste(directory, "Bar_gene_property_developmental_analysis_female_male.png", sep="")
ggsave(plot, width = 24, height = 10, unit = "cm", file = outfile)
