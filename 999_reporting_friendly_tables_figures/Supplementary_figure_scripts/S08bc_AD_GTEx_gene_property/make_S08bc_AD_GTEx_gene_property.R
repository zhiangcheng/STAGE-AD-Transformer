# Create supplementary fugure - Plots of gene-property analysis run in FUMA using sex-stratified GWAS sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Plot_gene_property_female_male/"

Magma_30_tissues_female = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"

Magma_30_tissues_male = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Males_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"

Magma_53_tissues_female = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"

Magma_53_tissues_male = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_Males_genomewidesig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"



### Packages ###
library(patchwork)
library(tidyverse)

### Load Data ###
Magma_30_tissues_female_df <- read.table(Magma_30_tissues_female, header = T, stringsAsFactors = F)

Magma_30_tissues_male_df <- read.table(Magma_30_tissues_male, header = T, stringsAsFactors = F)

Magma_53_tissues_female_df <- read.table(Magma_53_tissues_female, header = T, stringsAsFactors = F)

Magma_53_tissues_male_df <- read.table(Magma_53_tissues_male, header = T, stringsAsFactors = F)


### Dataframe: 30 tissues ### 

# Combine the data from both GWAS' into a single dataframe
combined_30_df <- bind_rows(
  mutate(Magma_30_tissues_female_df, Dataset = "female"),
  mutate(Magma_30_tissues_male_df, Dataset = "male")
) %>% 
  mutate(VARIABLE = str_replace(VARIABLE, "_", " "))


# Reorder the VARIABLE factor based on the P values for the female dataset
female_order <- combined_30_df %>%
  filter(Dataset == "female") %>%
  arrange(P) %>% 
  pull(VARIABLE)

# Apply the reordered factor levels to the whole dataset
combined_30_df <- combined_30_df %>%
  mutate(VARIABLE = factor(VARIABLE, levels = female_order))


### Plot ###

plot_30_tissues <- ggplot(data = combined_30_df, aes(x = VARIABLE, y = -log10(P), fill = Dataset)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  geom_hline(yintercept = -log10(0.05/30),
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

plot_30_tissues



### Dataframe: 53 tissues ### 

# Combine the data from both GWAS' into a single dataframe
combined_53_df <- bind_rows(
  mutate(Magma_53_tissues_female_df, Dataset = "female"),
  mutate(Magma_53_tissues_male_df, Dataset = "male")
) %>% 
  mutate(FULL_NAME = str_replace_all(FULL_NAME, "_", " "),
         FULL_NAME = str_wrap(FULL_NAME, width = 30))


# Reorder the FULL_NAME factor based on the P values for the female dataset
female_order <- combined_53_df %>%
  filter(Dataset == "female") %>%
  arrange(P) %>% 
  pull(FULL_NAME)

# Apply the reordered factor levels to the whole dataset
combined_53_df <- combined_53_df %>%
  mutate(FULL_NAME = factor(FULL_NAME, levels = female_order))


### Plot ###

plot_53_tissues <- ggplot(data = combined_53_df, aes(x = FULL_NAME, y = -log10(P), fill = Dataset)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d("Sex", begin = 1, end = 0,
                       labels = c("Females", "Males")) +
  geom_hline(yintercept = -log10(0.05/53),
             colour = "grey30",
             linetype = 'dotted') +
  labs(x = "Tissue",
       y = expression(-log[10](p-value))) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        legend.position = "none",
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))

plot_53_tissues



### Panelled Figure ###
figure <- plot_30_tissues /
  plot_53_tissues +
  plot_annotation(tag_levels = 'A') 

figure


outfile <- paste(directory, "Bar_gene_property_analysis_female_male.png", sep="")
ggsave(figure, width = 28, height = 20, unit = "cm", file = outfile)



