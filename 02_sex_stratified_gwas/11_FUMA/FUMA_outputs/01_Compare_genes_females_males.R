# Compare genes annotated from female and male genome wide significant SNPs

directory = "/path/07_FUMA/Downloaded_results/"

male = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Males_genomewidesig_positionalmap10kb_MHCexcl/Males_sig_SNPs_annotations_nearest_twoplusmethods.txt"

female="/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/Females_sig_SNPs_annotations_nearest_twoplusmethods.txt"

### LOAD PACKAGES ###

library(tidyverse)

### LOAD DATA ###

male_genes_df <- read.table(male, header = T, stringsAsFactors = F, sep = "\t")

female_genes_df <- read.table(female, header = T, stringsAsFactors = F, sep = "\t")

### COMPARE GENES FROM FEMALES AND MALES ###

genes_male <- male_genes_df %>% 
  rename(Genes_male = Genes..annotated.with...1.method.) %>% 
  select(Genes_male) %>%
  separate_rows(Genes_male, sep = ", ") %>% 
  distinct(Genes_male) %>% 
  mutate(Genes = Genes_male)

genes_female <- female_genes_df %>% 
  rename(Genes_female = Genes..annotated.with...1.method.) %>% 
  select(Genes_female) %>%
  mutate(Genes_female = na_if(Genes_female, "")) %>% 
  drop_na(Genes_female) %>% 
  separate_rows(Genes_female, sep = ", ") %>% 
  distinct(Genes_female) %>% 
  mutate(Genes = Genes_female)

genes_all <- genes_male %>% 
  full_join(genes_female, by = "Genes") %>% 
  select(Genes, Genes_female, Genes_male) %>% 
  arrange(Genes) %>% 
  mutate(
    Females = ifelse(is.na(Genes_female), "×", "✓"),
    Males = ifelse(is.na(Genes_male), "×", "✓")
  ) %>% 
  select(Genes, Females, Males)

# Save
outfile <- paste(directory, "Genes_annotated_with_more_than_one_method_females_males.txt", sep="")
write.table(genes_all, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

  
  
