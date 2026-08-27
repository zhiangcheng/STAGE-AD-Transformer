# Create tables of gene annotations for female GWAS

directory = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/"

snps_female = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/snps.txt"

genes_female = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_genomewidesig_positionalmap10kb_MHCexcl/genes.txt"


femaleGWAS="/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

maleGWAS="/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"

GxS="/path/03_Metal/GxS/fulldc/Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID.txt"


### LOAD PACKAGES ###

library(tidyverse)

### LOAD DATA ###

snps_female_df <- read.table(snps_female, header = T, stringsAsFactors = F)

genes_female_df <- read.table(genes_female, header = T, stringsAsFactors = F)

females_GWAS_df <- read.table(femaleGWAS, header = T, stringsAsFactors = F)

males_GWAS_df <- read.table(maleGWAS, header = T, stringsAsFactors = F)

GxS_df <- read.table(GxS, header = T, stringsAsFactors = F)



### CREATE TABLE OF GENOME-WIDE SIG SNPS FROM FEMALE GWAS, FEMALE GWAS RESULTS, MALE GWAS RESULTS, GxS RESULTS, ALL MAPPED GENES ###
# Independent lead SNP and each of the  sig SNPs in LD with it
snps_female_LD <- snps_female_df %>%
  group_by(IndSigSNP) %>%
  summarize(sig_SNPs_in_LD = toString(na.omit(rsID[rsID != IndSigSNP]))) %>%
  ungroup() %>%
  rename(IndSigSNPs = IndSigSNP)

# Mapped genes to each independent sig SNP (and the sig SNPs in LD with it) - each mapping type
genes_female_eqtl <- genes_female_df %>%
  mutate(eqtl_gene = case_when(
    eqtlMapSNPs > 0 ~ paste0(symbol, "(", eqtlMapts, ")")
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(eqtl_gene_tissue = toString(unique(na.omit(eqtl_gene)))) %>%
  ungroup()

genes_female_ci <- genes_female_df %>%
  mutate(ci_gene = case_when(
    ciMap == "Yes" ~ paste0(symbol, "(", ciMapts, ")")
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(ci_gene_tissue = toString(unique(na.omit(ci_gene)))) %>%
  ungroup()

genes_female_pos <- genes_female_df %>%
  mutate(pos_gene = case_when(
    posMapSNPs > 0 ~ paste0(symbol)
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(pos_gene = toString(unique(na.omit(pos_gene)))) %>%
  ungroup()


# Now merge into one dataset
genes_female_all <- snps_female_LD %>%
  full_join(genes_female_pos) %>%
  full_join(genes_female_eqtl) %>%
  full_join(genes_female_ci)


# Add females GWAS info
females_GWAS_df2 <- females_GWAS_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

snps_genes_female <- genes_female_all %>%
  left_join(females_GWAS_df2, by = c("IndSigSNPs" = "rsID_build37"))

snps_genes_female <- snps_genes_female %>%
  select(IndSigSNPs, CHR, BP, Allele1, Allele2, Effect, StdErr, P, Freq1, sig_SNPs_in_LD, pos_gene, eqtl_gene_tissue, ci_gene_tissue) %>%
  rename(rsID = IndSigSNPs,
         Chr = CHR,
         Pos = BP,
         "Effect Allele" = Allele1,
         "Other Allele" = Allele2,
         "Beta (Females)" = Effect,
         "SE (Females)" = StdErr,
         "P (Females)" = P,
         "Freq Effect Allele (Females)" = Freq1,
         "Significant SNPs in LD" = sig_SNPs_in_LD,
         "Gene: positional" = pos_gene,
         "Gene: eQTL (tissue)" = eqtl_gene_tissue,
         "Gene: chromatin interaction (tissue)" = ci_gene_tissue)


# Add male GWAS sumstats info
males_GWAS_df2 <- males_GWAS_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

snps_genes_female <- snps_genes_female %>%
  left_join(males_GWAS_df2, by = c("rsID" = "rsID_build37", "Effect Allele" = "Allele1", "Other Allele" = "Allele2"))

snps_genes_female <- snps_genes_female %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (Females)", "SE (Females)", "P (Females)", "Freq Effect Allele (Females)",
         Effect, StdErr, P, Freq1,
         "Significant SNPs in LD", "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)") %>%
  rename("Beta (Males)" = Effect,
         "SE (Males)" = StdErr,
         "P (Males)" = P,
         "Freq Effect Allele (Males)" = Freq1)


# Add GxS sumstats info
GxS_df2 <- GxS_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

snps_genes_female <- snps_genes_female %>%
  left_join(GxS_df2, by = c("rsID" = "rsID_build37", "Effect Allele" = "Allele1", "Other Allele" = "Allele2"))

snps_genes_female <- snps_genes_female %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (Females)", "SE (Females)", "P (Females)", "Freq Effect Allele (Females)",
         "Beta (Males)", "SE (Males)", "P (Males)", "Freq Effect Allele (Males)",
         Effect, StdErr, P,
         "Significant SNPs in LD", "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)") %>%
  rename("Beta (GxS)" = Effect,
         "SE (GxS)" = StdErr,
         "P (GxS)" = P)

# Save
outfile <- paste(directory, "Females_sig_SNPs_Femalesinfo_Malesinfo_GxSinfo_annotations.txt", sep="")
write.table(snps_genes_female, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



### CREATE TABLE OF GENOME-WIDE SIG FEMALE SNPs AND GENES MAPPED WITH > 1 METHOD ###
# For each independent SNP (and sig SNPs in LD with it) find genes mapped using more than 1 method
genes_2plus_methods_female <- genes_female_df %>%
  mutate(gene_2plus_methods = case_when(
    eqtlMapSNPs > 0 & ciMap == "Yes" ~ symbol,
    eqtlMapSNPs > 0 & posMapSNPs > 0 ~ symbol,
    ciMap == "Yes" & posMapSNPs > 0 ~ symbol
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(gene_2plus_methods = toString(unique(na.omit(gene_2plus_methods)))) %>%
  ungroup() %>%
  rename(rsID = IndSigSNPs,
         "Genes (annotated with > 1 method)" = gene_2plus_methods)

# Identify nearest gene (and its distance and function) for lead SNPs and sig SNPs in LD with lead SNPs
genes_nearest_snps_LD_female <- snps_female_df %>%
  mutate(nearest_gene_sig_SNPs_in_LD_dist_func = paste0(nearestGene, "(", dist, ", ", func, ")")) %>%
  group_by(IndSigSNP) %>%
  summarize(nearest_gene_sig_SNPs_in_LD_dist_func = toString(unique(na.omit(nearest_gene_sig_SNPs_in_LD_dist_func[rsID != IndSigSNP])),)) %>%
  ungroup() %>%
  rename(rsID = IndSigSNP,
         "Gene nearest significant SNPs in LD (distance, function)" = nearest_gene_sig_SNPs_in_LD_dist_func)

genes_nearest_lead_snps_female <- snps_female_df %>%
  filter(IndSigSNP == rsID) %>%
  mutate(nearest_gene_lead_SNP = paste0(nearestGene, "(", dist, ", ", func, ")")) %>%
  select(rsID, chr, pos, effect_allele, non_effect_allele, nearest_gene_lead_SNP) %>%
  rename(Chr = chr,
         Pos = pos,
         "Effect Allele" = effect_allele,
         "Other Allele" = non_effect_allele,
         "Gene nearest lead SNP (distance, function)" = nearest_gene_lead_SNP)

# Combine dataframes
genes_female_combined <- genes_nearest_lead_snps_female %>%
  full_join(genes_nearest_snps_LD_female) %>%
  full_join(genes_2plus_methods_female)

# Save
outfile <- paste(directory, "Females_sig_SNPs_annotations_nearest_twoplusmethods.txt", sep="")
write.table(genes_female_combined, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



### CREATE TABLE OF GENOME-WIDE SIG FEMALE SNPs, GWAS RESULTS, GENES MAPPED WITH > 1 METHOD, POS, eQTL, CI MAPPING ###

# Combine above two tables
female_large_table <- snps_genes_female %>%
  full_join(genes_female_combined) %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (Females)", "SE (Females)", "P (Females)", "Freq Effect Allele (Females)",
         "Beta (Males)", "SE (Males)", "P (Males)", "Freq Effect Allele (Males)",
         "Beta (GxS)", "SE (GxS)", "P (GxS)",
         "Significant SNPs in LD",
         "Gene nearest lead SNP (distance, function)", "Gene nearest significant SNPs in LD (distance, function)",
         "Genes (annotated with > 1 method)",
         "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)")

# Save
outfile <- paste(directory, "Females_sig_SNPs_GWAS_results_annotations_nearest_twoplusmethods_pos_eqtl_ci.txt", sep="")
write.table(female_large_table, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

