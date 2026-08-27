# Create supplementary figure - Plots of gene-property analysis run in FUMA using GxS sumstats

directory = "/path/10_Publication_figures/Supplementary_figures/Plot_gene_property_GxS_fulldc/"

Magma_30_tissues_GxS = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"

Magma_53_tissues_GxS = "/path/07_FUMA/Downloaded_results/AD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"




### Packages ###
library(patchwork)
library(tidyverse)

### Load Data ###
Magma_30_tissues_GxS_df <- read.table(Magma_30_tissues_GxS, header = T, stringsAsFactors = F)

Magma_53_tissues_GxS_df <- read.table(Magma_53_tissues_GxS, header = T, stringsAsFactors = F)




### Dataframe: 30 tissues ### 

Magma_30_tissues_GxS_df <- Magma_30_tissues_GxS_df %>% 
  mutate(VARIABLE = str_replace_all(VARIABLE, "_", " "))


# Reorder the VARIABLE factor based on the P values 
order <- Magma_30_tissues_GxS_df %>%
  arrange(P) %>% 
  pull(VARIABLE)

# Apply the reordered factor levels to the whole dataset
Magma_30_tissues_GxS_df <- Magma_30_tissues_GxS_df %>%
  mutate(VARIABLE = factor(VARIABLE, levels = order))


### Plot ###

plot_30_tissues <- ggplot(data = Magma_30_tissues_GxS_df, aes(x = VARIABLE, y = -log10(P))) +
  geom_col(position = position_dodge(),
           fill = "#20A387FF") +
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

### Dataframe: 30 tissues ### 

Magma_30_tissues_GxS_df <- Magma_30_tissues_GxS_df %>% 
  mutate(VARIABLE = str_replace_all(VARIABLE, "_", " "))


# Reorder the VARIABLE factor based on the P values 
order <- Magma_30_tissues_GxS_df %>%
  arrange(P) %>% 
  pull(VARIABLE)

# Apply the reordered factor levels to the whole dataset
Magma_30_tissues_GxS_df <- Magma_30_tissues_GxS_df %>%
  mutate(VARIABLE = factor(VARIABLE, levels = order))


### Plot ###

plot_30_tissues <- ggplot(data = Magma_30_tissues_GxS_df, aes(x = VARIABLE, y = -log10(P))) +
  geom_col(position = position_dodge(),
           fill = "#20A387FF") +
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

Magma_53_tissues_GxS_df <- Magma_53_tissues_GxS_df %>% 
  mutate(FULL_NAME = str_replace_all(FULL_NAME, "_", " "),
         FULL_NAME = str_wrap(FULL_NAME, width = 30))


# Reorder the FULL_NAME factor based on the P values 
order <- Magma_53_tissues_GxS_df %>%
  arrange(P) %>% 
  pull(FULL_NAME)

# Apply the reordered factor levels to the whole dataset
Magma_53_tissues_GxS_df <- Magma_53_tissues_GxS_df %>%
  mutate(FULL_NAME = factor(FULL_NAME, levels = order))


### Plot ###

plot_53_tissues <- ggplot(data = Magma_53_tissues_GxS_df, aes(x = FULL_NAME, y = -log10(P))) +
  geom_col(position = position_dodge(),
           fill = "#20A387FF") +
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


outfile <- paste(directory, "Bar_gene_property_analysis_GxS_fulldc.png", sep="")
ggsave(figure, width = 28, height = 20, unit = "cm", file = outfile)



