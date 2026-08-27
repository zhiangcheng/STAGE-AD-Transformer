# Compare male and female polygenicity (pi) from meta-analysed sumstats and after combining all per cohort pi estimates


directory = "/path/12_SBayesS_h2/"


MA_females_mcmc <- paste0(directory, "Metaanalysis_MDD_female_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

MA_males_mcmc <- paste0(directory, "Metaanalysis_MDD_male_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")


AGDS_females_mcmc <- paste0(directory, "AGDS_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
AGDS_males_mcmc <- paste0(directory, "AGDS_MDD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

AllOfUs_females_mcmc <- paste0(directory, "AllOfUs_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
AllOfUs_males_mcmc <- paste0(directory, "AllOfUs_MDD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

Bionic_females_mcmc <- paste0(directory, "Bionic_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
Bionic_males_mcmc <- paste0(directory, "Bionic_MDD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

Blokland21_females_mcmc <- paste0(directory, "Blokland21_MDD_female_sumstats_SBayesS_DENTISTremoved_redo.mcmcsamples/CoreParameters.mcmcsamples.txt")
# update when run done properly
Blokland21_males_mcmc <- paste0(directory, "Blokland21_MDD_male_sumstats_trimmedN_30perc_SBayesS.mcmcsamples/CoreParameters.mcmcsamples.txt")

GLAD_females_mcmc <- paste0(directory, "GLAD_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
# update when run done properly
GLAD_males_mcmc <- paste0(directory, "GLAD_MDD_male_sumstats_SBayesS_h006_pi001.mcmcsamples/CoreParameters.mcmcsamples.txt")

UKB_females_mcmc <- paste0(directory, "UKB_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
UKB_males_mcmc <- paste0(directory, "UKB_MDD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")



### Packages ###

library(metafor)      # meta-analysis
library(coda)         # HPDI
library(tidyverse)



### Load Data ###

# Load all mcmc reports

MA_females_mcmc_df <- read.table(MA_females_mcmc, header = TRUE, stringsAsFactors = FALSE)
MA_males_mcmc_df <- read.table(MA_males_mcmc, header = TRUE, stringsAsFactors = FALSE)

MA_females_mcmc_df <- MA_females_mcmc_df %>%
  rename_with(~ paste0(., "_female"))

MA_males_mcmc_df <- MA_males_mcmc_df %>%
  rename_with(~ paste0(., "_male"))

MA_all_mcmc <- bind_cols(MA_females_mcmc_df, MA_males_mcmc_df)

# Remove burn-in of 5,000 (as thinning of 10 used = 500 rows)
MA_all_mcmc <- MA_all_mcmc %>% 
  slice(501:n())


### Calculate mean Pi +- HPDI: GWAS M-A sumstats ###

# Female
posterior_mean_pi_female <- MA_all_mcmc %>%
  summarize(mean = mean(Pi_female)) %>%
  pull(mean)
hpd_interval_pi_female <- HPDinterval(as.mcmc(MA_all_mcmc %>% pull(Pi_female)), prob = 0.95)[1, ]

# Male
posterior_mean_pi_male <- MA_all_mcmc %>%
  summarize(mean = mean(Pi_male)) %>%
  pull(mean)
hpd_interval_pi_male <- HPDinterval(as.mcmc(MA_all_mcmc %>% pull(Pi_male)), prob = 0.95)[1, ]


# Create data frame for values
pi_MA <- data.frame(
  mean_pi = c(posterior_mean_pi_female, posterior_mean_pi_male),
  hpdi_lower = c(hpd_interval_pi_female[1], hpd_interval_pi_male[1]),
  hpdi_upper = c(hpd_interval_pi_female[2], hpd_interval_pi_male[2]),
  sex = c("female", "male"),
  Type = c("MA", "MA")
)

### Male vs Female: GWAS M-A sumstats ###
# posterioir prob that pi female is > pi male

# Count the frequency in the MCMC sample that pi_female > pi_male = posterior probability that female pi is larger than male pi
post_prob_M_F <- MA_all_mcmc %>%
  mutate(pi_f_minus_m = Pi_female - Pi_male) %>%
  summarize(Perc = (sum(pi_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_Pi_GWAS_MA_sumstats.txt"))
cat("Posterior probability that Pi_female > Pi_male:", post_prob_M_F, "\n")
sink()



### LOAD DATA ALL COHORTS ###

read_and_rename_mcmc <- function(file_path, suffix) {
  read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    slice(501:n()) %>%
    rename_with(~ paste0(., "_", suffix))
}

# Process each cohort
AGDS_all_mcmc <- bind_cols(
  read_and_rename_mcmc(AGDS_females_mcmc, "female"),
  read_and_rename_mcmc(AGDS_males_mcmc, "male")
)

AllOfUs_all_mcmc <- bind_cols(
  read_and_rename_mcmc(AllOfUs_females_mcmc, "female"),
  read_and_rename_mcmc(AllOfUs_males_mcmc, "male")
)

Bionic_all_mcmc <- bind_cols(
  read_and_rename_mcmc(Bionic_females_mcmc, "female"),
  read_and_rename_mcmc(Bionic_males_mcmc, "male")
)

Blokland21_all_mcmc <- bind_cols(
  read_and_rename_mcmc(Blokland21_females_mcmc, "female"),
  read_and_rename_mcmc(Blokland21_males_mcmc, "male")
)

GLAD_all_mcmc <- bind_cols(
  read_and_rename_mcmc(GLAD_females_mcmc, "female"),
  read_and_rename_mcmc(GLAD_males_mcmc, "male")
)

UKB_all_mcmc <- bind_cols(
  read_and_rename_mcmc(UKB_females_mcmc, "female"),
  read_and_rename_mcmc(UKB_males_mcmc, "male")
)



# Create one df with data from all cohorts
Cohorts_all_mcmc <- bind_rows(
  AGDS_all_mcmc,
  AllOfUs_all_mcmc,
  Bionic_all_mcmc,
  Blokland21_all_mcmc,
  GLAD_all_mcmc,
  UKB_all_mcmc
)


### Calculate mean Pi +- HPDI all cohorts combined ###

# Female
posterior_mean_pi_female <- Cohorts_all_mcmc %>%
  summarize(mean = mean(Pi_female)) %>%
  pull(mean)
hpd_interval_pi_female <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(Pi_female)), prob = 0.95)[1, ]

# Male
posterior_mean_pi_male <- Cohorts_all_mcmc %>%
  summarize(mean = mean(Pi_male)) %>%
  pull(mean)
hpd_interval_pi_male <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(Pi_male)), prob = 0.95)[1, ]


# Create data frame for values
pi_cohorts <- data.frame(
  mean_pi = c(posterior_mean_pi_female, posterior_mean_pi_male),
  hpdi_lower = c(hpd_interval_pi_female[1], hpd_interval_pi_male[1]),
  hpdi_upper = c(hpd_interval_pi_female[2], hpd_interval_pi_male[2]),
  sex = c("female", "male"),
  Type = c("cohorts", "cohorts")
)


### Male vs Female: Using all cohorts ###
# posterior prob that pi female is > pi male

# Count the frequency in the MCMC sample that pi_female > pi_male = posterior probability that female pi is larger than male pi
post_prob_M_F <- Cohorts_all_mcmc %>%
  mutate(pi_f_minus_m = Pi_female - Pi_male) %>%
  summarize(Perc = (sum(pi_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_pi_Cohorts.txt"))
cat("Posterior probability that Pi_female > Pi_male (all cohorts combined):", post_prob_M_F, "\n")
sink()





### Figure Pi female and male, MA and all cohorts ###

all_data_pi <- bind_rows(pi_MA, pi_cohorts)
  
# Save the combined data frame
save(all_data_pi, file = paste0(directory, "Pi_data_MA_cohorts.RData"))
write.table(all_data_pi, file = paste0(directory, "Pi_data_MA_cohorts.txt"), col.names = T, row.names = F, sep = "\t", quote = F)


pi_plot <- ggplot(all_data_pi, aes(x = mean_pi, y = Type, colour = sex)) +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(xmin = hpdi_lower, xmax = hpdi_upper),
                 position = position_dodge(width = -0.8)) +
  scale_y_discrete("",
                   breaks = c("cohorts",
                              "MA"),
                   limits = c("cohorts",
                              "MA"),
                   labels = c("Each cohort",
                              "Meta-analysed sumstats")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  scale_x_continuous("Polygenicity (\u03c0)",
                     limits = c(0,0.031)) +
  annotate("text", x = 0.028, y = "MA", label = "100%", size = 3) +
  annotate("text", x = 0.028, y = "cohorts", label = "79%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "right")

pi_plot

ggsave(pi_plot, file = paste0(directory, "forest_plot_Pi_MA_cohorts.png"), width = 15, height = 10, unit = "cm")













