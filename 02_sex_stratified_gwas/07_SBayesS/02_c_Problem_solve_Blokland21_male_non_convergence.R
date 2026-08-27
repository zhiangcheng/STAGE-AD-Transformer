# Blokland male sumstats having problems converging in SBayesS even when having removed problem SNPs as identified by DENTIST
# So look into SNPs and remove those with much lower N

directory = "/path/12_SBayesS_h2/"

GWAS_Blokland21_male = paste0(directory, "Blokland21_MDD_male_sumstats_formatted_SBayes.ma", sep = "")



### Load Packages ###

library(tidyverse)


### Load Data ###

GWAS_Blokland21_male_df <- read.table(GWAS_Blokland21_male, header = T)



### Check per-SNP sample size ###

hist(GWAS_Blokland21_male_df$N)

GWAS_Blokland21_male_df %>%
  summarise(median = median(N),
            SD = sd(N),
            median_minus_3SD = median(N) - (3 * sd(N)),
            median_minus_2SD = median(N) - (2 * sd(N)),
            median_minus_1SD = median(N) - (1 * sd(N)),
            lowest_10perc = quantile(N, 0.1),
            lowest_20perc = quantile(N, 0.2),
            lowest_30perc = quantile(N, 0.3),
            lowest_40perc = quantile(N, 0.4)
  )


# Remove SNPs with N < 3 SD below the median

GWAS_Blokland21_male_trimmed_3SD_df <- GWAS_Blokland21_male_df %>%
  filter(N > 4908)

hist(GWAS_Blokland21_male_trimmed_3SD_df$N)


write.table(GWAS_Blokland21_male_trimmed_3SD_df, paste0(directory, "Blokland21_MDD_male_sumstats_formatted_trimmedN_3SD_SBayes.ma"), quote = F, col.names = T, row.names = F, sep = "\t")

# Remove SNPs in lowest 10% quantile

GWAS_Blokland21_male_trimmed_10perc_df <- GWAS_Blokland21_male_df %>%
  filter(N > 10684)

hist(GWAS_Blokland21_male_trimmed_10perc_df$N)


write.table(GWAS_Blokland21_male_trimmed_10perc_df, paste0(directory, "Blokland21_MDD_male_sumstats_formatted_trimmedN_10perc_SBayes.ma"), quote = F, col.names = T, row.names = F, sep = "\t")

# Remove SNPs in lowest 20% quantile

GWAS_Blokland21_male_trimmed_20perc_df <- GWAS_Blokland21_male_df %>%
  filter(N > 13858)

hist(GWAS_Blokland21_male_trimmed_20perc_df$N)


write.table(GWAS_Blokland21_male_trimmed_20perc_df, paste0(directory, "Blokland21_MDD_male_sumstats_formatted_trimmedN_20perc_SBayes.ma"), quote = F, col.names = T, row.names = F, sep = "\t")


# Remove SNPs in lowest 30% quantile

GWAS_Blokland21_male_trimmed_30perc_df <- GWAS_Blokland21_male_df %>%
  filter(N > 15871)

hist(GWAS_Blokland21_male_trimmed_30perc_df$N)


write.table(GWAS_Blokland21_male_trimmed_30perc_df, paste0(directory, "Blokland21_MDD_male_sumstats_formatted_trimmedN_30perc_SBayes.ma"), quote = F, col.names = T, row.names = F, sep = "\t")

# Remove SNPs in lowest 40% quantile

GWAS_Blokland21_male_trimmed_40perc_df <- GWAS_Blokland21_male_df %>%
  filter(N > 24372)

hist(GWAS_Blokland21_male_trimmed_40perc_df$N)


write.table(GWAS_Blokland21_male_trimmed_40perc_df, paste0(directory, "Blokland21_MDD_male_sumstats_formatted_trimmedN_40perc_SBayes.ma"), quote = F, col.names = T, row.names = F, sep = "\t")

