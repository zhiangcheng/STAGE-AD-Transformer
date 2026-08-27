# Check GWAS Catalog matches for top SNPs (and SNPs in LD) (FUMA output) and see if any novel or all previously shown to be associated with MDD.
# Use FUMA run that included all SNPs in LD with genome-wide significant SNPs (not just sig SNPs in LD) as well as SNPs in LD that aren't in my sumstats but are present in 1000Genomes

directory = "/path/07_FUMA/Downloaded_results/"

gwascatalog_female = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_allSNP_positional10kb_MHCincl/gwascatalog.txt"

snps_female = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Females_allSNP_positional10kb_MHCincl/snps.txt"

gwascatalog_male = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Males_allSNP_positional10kb_MHCincl/gwascatalog.txt"

snps_male = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_Males_allSNP_positional10kb_MHCincl/snps.txt"

gwascatalog_gxs = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_allSNP_positional10kb_MHCincl/gwascatalog.txt"

snps_gxs = "/path/07_FUMA/Downloaded_results/MDD_Metaanalysis_Freeze3_GxS_fulldc_redonetransformation_allSNP_positional10kb_MHCincl/snps.txt"

AdamsGWAS = "/path/Summary_statistics_published/Depression_Adams_2024/including_23andMe/Adams_Full_EUR_2025_FromIMB.txt"

Adams_clumping = "/path/09_Effect_sizes_plots/Adams_SNPs/clumped_Adams_SNPs_all_results.txt"



### LOAD DATA ###

# had problems reading in so use read.delim and works ok
gwascatalog_female_df <- read.delim(gwascatalog_female, header = T, stringsAsFactors = F, sep = "\t")

snps_female_df <- read.table(snps_female, header = T, stringsAsFactors = F)

gwascatalog_male_df <- read.delim(gwascatalog_male, header = T, stringsAsFactors = F, sep = "\t")

snps_male_df <- read.table(snps_male, header = T, stringsAsFactors = F)

gwascatalog_gxs_df <- read.delim(gwascatalog_gxs, header = T, stringsAsFactors = F, sep = "\t")

snps_gxs_df <- read.table(snps_gxs, header = T, stringsAsFactors = F)

AdamsGWAS_df <- read.table(AdamsGWAS, header = T, stringsAsFactors = F, comment.char = "")

Adams_clumping_df <- read.table(Adams_clumping, header = T, stringsAsFactors = F)



### FORMAT ADAMS GWAS ###
# latest GWAS of depression by Adams isn't yet in GWASCatalog so check if any of our hits are in there
# Create list of all genome-wide sig SNPs in Adams + any SNPs in LD with these genome-wide sig SNPs
Adams_sig_and_LD <- Adams_clumping_df %>% 
  mutate(sig_SNPs_and_LD = paste(SNP, SP2, sep = ","),
         sig_SNPs_and_LD = str_replace_all(sig_SNPs_and_LD, "\\(1\\)", "")) %>%  
  separate_longer_delim(sig_SNPs_and_LD, delim = ",") %>%    
  distinct(sig_SNPs_and_LD) %>%                                     
  select(sig_SNPs_and_LD)                                           

## FEMALES
# If any of the female MDD lead SNPs, or SNPS in LD, are also in Adams lead SNPs, or SNPs in LD then add trait "Adams Major Depressive Disorder"
matches_female <- snps_female_df %>%
  filter(rsID %in% Adams_sig_and_LD$sig_SNPs_and_LD) %>%
  mutate(Trait = "Adams Major Depressive Disorder") %>%
  select(GenomicLocus, IndSigSNP, chr, pos, rsID, Trait) %>%
  rename(bp = pos,
         snp = rsID)
# Now add these SNPs to gwascatalog df (will ensure capturing SNPs in Adams that weren't in gwascatalog df because the snp had no gwascatalog match)
gwascatalog_female_df2 <- gwascatalog_female_df %>%
  full_join(matches_female)


## MALES
# Add Adams sig SNPs as trait to any of the lead SNPs or SNPs in LD from my GWAS
matches_male <- snps_male_df %>%
  filter(rsID %in% Adams_sig_and_LD$sig_SNPs_and_LD) %>%
  mutate(Trait = "Adams Major Depressive Disorder") %>%
  select(GenomicLocus, IndSigSNP, chr, pos, rsID, Trait) %>%
  rename(bp = pos,
         snp = rsID)
# Now add these SNPs to gwascatalog df (will ensure capturing SNPs in Adams that weren't in gwascatalog df because the snp had no gwascatalog match)
gwascatalog_male_df2 <- gwascatalog_male_df %>%
  full_join(matches_male)



## GxS
# Add Adams sig SNPs as trait to any of the lead SNPs or SNPs in LD from my GWAS
matches_gxs <- snps_gxs_df %>%
  filter(rsID %in% Adams_sig_and_LD$sig_SNPs_and_LD) %>%
  mutate(Trait = "Adams Major Depressive Disorder") %>%
  select(GenomicLocus, IndSigSNP, chr, pos, rsID, Trait) %>%
  rename(bp = pos,
         snp = rsID)
# Now add these SNPs to gwascatalog df (will ensure capturing SNPs in Adams that weren't in gwascatalog df because the snp had no gwascatalog match)
gwascatalog_gxs_df2 <- gwascatalog_gxs_df %>%
  full_join(matches_gxs)



### FEMALES: LOOK AT GWASCATALOG HITS FOR EACH SIG SNP (AND SNPs IN LD) ###

# What traits have just the independent significant SNPs been mapped to?
indsigsnps <- gwascatalog_female_df2 %>%
  distinct(IndSigSNP) %>%
  pull(IndSigSNP)

indsigsnps_df <- gwascatalog_female_df2 %>%
  distinct(IndSigSNP)

traits_indsigsnps <- gwascatalog_female_df2 %>%
  filter(snp %in% indsigsnps) %>%
  select(snp, Trait) %>%
  arrange(Trait) %>%
  group_by(snp) %>%
  summarize(Trait = toString(unique(na.omit(Trait))), .groups = "drop") %>%
  right_join(indsigsnps_df, by = c("snp" = "IndSigSNP")) %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))


outfile <- paste(directory, "Females_ind_sig_snps_GWASCatalog_traits.txt", sep="")
write.table(traits_indsigsnps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# What traits have the independent significant SNPs and all the SNPs in LD with them been mapped to?
traits_snps <- gwascatalog_female_df2 %>%
  select(IndSigSNP, Trait) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait)))) %>%
  rename("Independent signficant SNP" = IndSigSNP)


outfile <- paste(directory, "Females_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits.txt", sep="")
write.table(traits_snps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# Does each independent snp (and snps in LD) map to a depression trait?
traits_snps_depr <- gwascatalog_female_df %>%
  select(IndSigSNP, Trait) %>%
  filter(str_detect(Trait, "(?i)depr")) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait)))) %>%
  right_join(indsigsnps_df, by = "IndSigSNP") %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))

# Pull out lead SNPs in which none of them or their SNPs in LD map to a depressive trait
snps_no_depr <- traits_snps_depr %>%
  filter(is.na(Trait)) %>%
  pull(IndSigSNP)

snps_no_depr_df <- traits_snps %>%
  filter(`Independent signficant SNP` %in% snps_no_depr)

outfile <- paste(directory, "Females_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits_no_depression.txt", sep="")
write.table(snps_no_depr_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)





### MALES: LOOK AT GWASCATALOG HITS FOR EACH SIG SNP (AND SNPs IN LD) ###

# What traits have just the independent significant SNPs been mapped to?
indsigsnps <- gwascatalog_male_df2 %>%
  distinct(IndSigSNP) %>%
  pull(IndSigSNP)

indsigsnps_df <- gwascatalog_male_df2 %>%
  distinct(IndSigSNP)

traits_indsigsnps <- gwascatalog_male_df2 %>%
  filter(snp %in% indsigsnps) %>%
  select(snp, Trait) %>%
  arrange(Trait) %>%
  group_by(snp) %>%
  summarize(Trait = toString(unique(na.omit(Trait))), .groups = "drop") %>%
  right_join(indsigsnps_df, by = c("snp" = "IndSigSNP")) %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))


outfile <- paste(directory, "Males_ind_sig_snps_GWASCatalog_traits.txt", sep="")
write.table(traits_indsigsnps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# What traits have the independent significant SNPs and all the SNPs in LD with them been mapped to?
traits_snps <- gwascatalog_male_df2 %>%
  select(IndSigSNP, Trait) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait)))) %>%
  rename("Independent signficant SNP" = IndSigSNP)


outfile <- paste(directory, "Males_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits.txt", sep="")
write.table(traits_snps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# Does each independent snp (and snps in LD) map to a depression trait?
traits_snps_depr <- gwascatalog_male_df2 %>%
  select(IndSigSNP, Trait) %>%
  filter(str_detect(Trait, "(?i)depr")) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait)))) %>%
  right_join(indsigsnps_df, by = "IndSigSNP") %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))

# Pull out lead SNPs in which none of them or their SNPs in LD map to a depressive trait
snps_no_depr <- traits_snps_depr %>%
  filter(is.na(Trait)) %>%
  pull(IndSigSNP)

snps_no_depr_df <- traits_snps %>%
  filter(`Independent signficant SNP` %in% snps_no_depr)

outfile <- paste(directory, "Males_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits_no_depression.txt", sep="")
write.table(snps_no_depr_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)





### GxS: LOOK AT GWASCATALOG HITS FOR EACH SIG SNP (AND SNPs IN LD) ###

# What traits have just the independent significant SNPs been mapped to? (nominally significant independent SNPs for GxS)
indsigsnps_df <- gwascatalog_gxs_df2 %>%
  distinct(IndSigSNP)

# Only pulls out 3 SNPs but there were 4 independent nominally significant SNPs from GxS analysis
# This means one has never been associated with a trait in GWAS catalog (itself or any SNPs in LD)
# Add this SNP to come up as NA (rs6080675)
new_row <- data.frame(IndSigSNP = "rs6080675")
indsigsnps_df <- indsigsnps_df %>%
  bind_rows(new_row)

indsigsnps <- indsigsnps_df %>%
  pull(IndSigSNP)


traits_indsigsnps <- gwascatalog_gxs_df2 %>%
  filter(snp %in% indsigsnps) %>%
  select(snp, Trait) %>%
  group_by(snp) %>%
  summarize(Trait = toString(unique(na.omit(Trait))), .groups = "drop") %>%
  right_join(indsigsnps_df, by = c("snp" = "IndSigSNP")) %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))

outfile <- paste(directory, "GxS_ind_sig_snps_GWASCatalog_traits.txt", sep="")
write.table(traits_indsigsnps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# What traits have the independent significant SNPs and all the SNPs in LD with them been mapped to?
traits_snps <- gwascatalog_gxs_df2 %>%
  select(IndSigSNP, Trait) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait))))%>%
  right_join(indsigsnps_df, by = "IndSigSNP") %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait)) %>%
  rename("Independent signficant SNP" = IndSigSNP)


outfile <- paste(directory, "GxS_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits.txt", sep="")
write.table(traits_snps, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



# Does each independent snp (and snps in LD) map to a depression trait?
traits_snps_depr <- gwascatalog_gxs_df2 %>%
  select(IndSigSNP, Trait) %>%
  filter(str_detect(Trait, "(?i)depr")) %>%
  arrange(Trait) %>%
  group_by(IndSigSNP) %>%
  distinct(Trait) %>%
  summarize(Trait = toString(unique(na.omit(Trait)))) %>%
  right_join(indsigsnps_df, by = "IndSigSNP") %>%
  mutate(Trait = ifelse(Trait == "", NA, Trait))

# Pull out lead SNPs in which none of them or their SNPs in LD map to a depressive trait
snps_no_depr <- traits_snps_depr %>%
  filter(is.na(Trait)) %>%
  pull(IndSigSNP)

snps_no_depr_df <- traits_snps %>%
  filter(`Independent signficant SNP` %in% snps_no_depr)

outfile <- paste(directory, "GxS_ind_sig_snps_and_all_snps_in_LD_GWASCatalog_traits_no_depression.txt", sep="")
write.table(snps_no_depr_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
