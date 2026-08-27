# Create a file that has all the info for the nominally significant SNPs from GxS analysis
# GxS info, Females GWAS info, Males GWAS info, nominally sig SNPs in LD with lead SNPs, positional, eQTL and chromatin interaction mappings


directory = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/"

snps = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/snps.txt"

genes = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_nominalsig_positionalmap10kb_MHCexcl/genes.txt"

femaleGWAS="/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

maleGWAS="/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"


### LOAD PACKAGES ###

library(tidyverse)

### LOAD DATA ###

snps_df <- read.table(snps, header = T, stringsAsFactors = F)

genes_df <- read.table(genes, header = T, stringsAsFactors = F)

femaleGWAS_df <- read.table(femaleGWAS, header = T, stringsAsFactors = F)

maleGWAS_df <- read.table(maleGWAS, header = T, stringsAsFactors = F)

### CREATE TABLE OF NOM SIG SNPS, GxS, FEMALE GWAS, MALE GWAS, ALL MAPPED GENES ###
# Independent lead SNP and each of the nominally sig SNPs in LD with it
snps_LD <- snps_df %>%
  group_by(IndSigSNP) %>%
  summarize(sig_SNPs_in_LD = toString(na.omit(rsID[rsID != IndSigSNP]))) %>%
  ungroup() %>%
  rename(IndSigSNPs = IndSigSNP)

# Mapped genes to each independent nominally sig SNP (and the nominally sig SNPs in LD with it) - each mapping type
genes_eqtl <- genes_df %>%
  mutate(eqtl_gene = case_when(
    eqtlMapSNPs > 0 ~ paste0(symbol, "(", eqtlMapts, ")")
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(eqtl_gene_tissue = toString(na.omit(eqtl_gene))) %>%
  ungroup()

genes_ci <- genes_df %>%
  mutate(ci_gene = case_when(
    ciMap == "Yes" ~ paste0(symbol, "(", ciMapts, ")")
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(ci_gene_tissue = toString(na.omit(ci_gene))) %>%
  ungroup()

genes_pos <- genes_df %>%
  mutate(pos_gene = case_when(
    posMapSNPs > 0 ~ paste0(symbol)
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(pos_gene = toString(na.omit(pos_gene))) %>%
  ungroup()


# Now merge into one dataset
genes_all <- snps_LD %>%
  full_join(genes_pos) %>%
  full_join(genes_eqtl) %>%
  full_join(genes_ci)


# Add GxS sumstats info
snps_genes <- genes_all %>%
  left_join(snps_df, by = c("IndSigSNPs" = "rsID")) %>%
  select(IndSigSNPs, chr, pos, effect_allele, non_effect_allele, beta, se, gwasP, sig_SNPs_in_LD, pos_gene, eqtl_gene_tissue, ci_gene_tissue) %>%
  rename(rsID = IndSigSNPs,
         Chr = chr,
         Pos = pos,
         "Effect Allele" = effect_allele,
         "Other Allele" = non_effect_allele,
         "Beta (GxS)" = beta,
         "SE (GxS)" = se,
         "P (GxS)" = gwasP,
         "Significant SNPs in LD" = sig_SNPs_in_LD,
         "Gene: positional" = pos_gene,
         "Gene: eQTL (tissue)" = eqtl_gene_tissue,
         "Gene: chromatin interaction (tissue)" = ci_gene_tissue)

# Add females GWAS info
femaleGWAS_df2 <- femaleGWAS_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

snps_genes <- snps_genes %>%
  left_join(femaleGWAS_df2, by = c("rsID" = "rsID_build37", "Effect Allele" = "Allele1", "Other Allele" = "Allele2"))

snps_genes <- snps_genes %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (GxS)", "SE (GxS)", "P (GxS)",
         Effect, StdErr, P, Freq1,
         "Significant SNPs in LD", "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)") %>%
  rename("Beta (Females)" = Effect,
         "SE (Females)" = StdErr,
         "P (Females)" = P,
         "Freq Effect Allele (Females)" = Freq1)



# Add male GWAS sumstats info
maleGWAS_df2 <- maleGWAS_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

snps_genes <- snps_genes %>%
  left_join(maleGWAS_df2, by = c("rsID" = "rsID_build37", "Effect Allele" = "Allele1", "Other Allele" = "Allele2"))

snps_genes <- snps_genes %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (GxS)", "SE (GxS)", "P (GxS)",
         "Beta (Females)", "SE (Females)", "P (Females)", "Freq Effect Allele (Females)",
         Effect, StdErr, P, Freq1,
         "Significant SNPs in LD", "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)") %>%
  rename("Beta (Males)" = Effect,
         "SE (Males)" = StdErr,
         "P (Males)" = P,
         "Freq Effect Allele (Males)" = Freq1)

# Save
outfile <- paste(directory, "GxS_nominally_sig_SNPs_GxSinfo_Femalesinfo_Malesinfo_annotations.txt", sep="")
write.table(snps_genes, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



### CREATE TABLE OF NOM SIG GxS SNPs AND GENES MAPPED WITH > 1 METHOD ###
# For each independent SNP (and sig SNPs in LD with it) find genes mapped using more than 1 method
genes_2plus_methods <- genes_df %>%
  mutate(gene_2plus_methods = case_when(
    eqtlMapSNPs > 0 & ciMap == "Yes" ~ symbol,
    eqtlMapSNPs > 0 & posMapSNPs > 0 ~ symbol,
    ciMap == "Yes" & posMapSNPs > 0 ~ symbol
  )) %>%
  group_by(IndSigSNPs) %>%
  summarize(gene_2plus_methods = toString(na.omit(gene_2plus_methods))) %>%
  ungroup() %>%
  rename(rsID = IndSigSNPs,
         "Genes (annotated with > 1 method)" = gene_2plus_methods)

# Identify nearest gene (and its distance and function) for lead SNPs and sig SNPs in LD with lead SNPs
genes_nearest_snps_LD <- snps_df %>%
  mutate(nearest_gene_sig_SNPs_in_LD_dist_func = paste0(nearestGene, "(", dist, ", ", func, ")")) %>%
  group_by(IndSigSNP) %>%
  summarize(nearest_gene_sig_SNPs_in_LD_dist_func = toString(na.omit(nearest_gene_sig_SNPs_in_LD_dist_func[rsID != IndSigSNP])),) %>%
  ungroup() %>%
  rename(rsID = IndSigSNP,
         "Gene nearest significant SNPs in LD (distance, function)" = nearest_gene_sig_SNPs_in_LD_dist_func)

genes_nearest_lead_snps <- snps_df %>%
  filter(IndSigSNP == rsID) %>%
  mutate(nearest_gene_lead_SNP = paste0(nearestGene, "(", dist, ", ", func, ")")) %>%
  select(rsID, chr, pos, effect_allele, non_effect_allele, nearest_gene_lead_SNP) %>%
  rename(Chr = chr,
         Pos = pos,
         "Effect Allele" = effect_allele,
         "Other Allele" = non_effect_allele,
         "Gene nearest lead SNP (distance, function)" = nearest_gene_lead_SNP)

# Combine dataframes
genes_combined <- genes_nearest_lead_snps %>%
  full_join(genes_nearest_snps_LD) %>%
  full_join(genes_2plus_methods)

# Save
outfile <- paste(directory, "GxS_nominally_sig_SNPs_annotations_nearest_twoplusmethods.txt", sep="")
write.table(genes_combined, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


### CREATE TABLE OF NOMINAL SIG GxS SNPs, GWAS RESULTS, GENES MAPPED WITH > 1 METHOD, POS, eQTL, CI MAPPING ###

# Combine above two tables
large_table <- snps_genes %>%
  full_join(genes_combined) %>%
  select(rsID, Chr, Pos, "Effect Allele", "Other Allele",
         "Beta (GxS)", "SE (GxS)", "P (GxS)",
         "Beta (Females)", "SE (Females)", "P (Females)", "Freq Effect Allele (Females)",
         "Beta (Males)", "SE (Males)", "P (Males)", "Freq Effect Allele (Males)",
         "Significant SNPs in LD",
         "Gene nearest lead SNP (distance, function)", "Gene nearest significant SNPs in LD (distance, function)",
         "Genes (annotated with > 1 method)",
         "Gene: positional", "Gene: eQTL (tissue)", "Gene: chromatin interaction (tissue)")

# Save
outfile <- paste(directory, "GxS_nominally_sig_SNPs_GWAS_results_annotations_nearest_twoplusmethods_pos_eqtl_ci.txt", sep="")
write.table(large_table, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
