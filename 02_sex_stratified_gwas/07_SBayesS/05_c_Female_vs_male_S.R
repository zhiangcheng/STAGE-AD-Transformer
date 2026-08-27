# Compare male and female selection factor (S) from meta-analysed sumstats and after combining all per cohort S estimates


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


### Calculate mean S +- HPDI: GWAS M-A sumstats ###

# Female
posterior_mean_S_female <- MA_all_mcmc %>%
  summarize(mean = mean(S_female)) %>%
  pull(mean)
hpd_interval_S_female <- HPDinterval(as.mcmc(MA_all_mcmc %>% pull(S_female)), prob = 0.95)[1, ]

# Male
posterior_mean_S_male <- MA_all_mcmc %>%
  summarize(mean = mean(S_male)) %>%
  pull(mean)
hpd_interval_S_male <- HPDinterval(as.mcmc(MA_all_mcmc %>% pull(S_male)), prob = 0.95)[1, ]


# Create data frame for values
S_MA <- data.frame(
  mean_S = c(posterior_mean_S_female, posterior_mean_S_male),
  hpdi_lower = c(hpd_interval_S_female[1], hpd_interval_S_male[1]),
  hpdi_upper = c(hpd_interval_S_female[2], hpd_interval_S_male[2]),
  sex = c("female", "male"),
  Type = c("MA", "MA")
)

### Male vs Female: GWAS M-A sumstats ###
# posterioir prob that S female is > S male

# Count the frequency in the MCMC sample that S_female > S_male = posterior probability that female S is larger than male S
post_prob_M_F <- MA_all_mcmc %>%
  mutate(S_f_minus_m = S_female - S_male) %>%
  summarize(Perc = (sum(S_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_S_GWAS_MA_sumstats.txt"))
cat("Posterior probability that S_female > S_male:", post_prob_M_F, "\n")
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


### Calculate mean S +- HPDI all cohorts combined ###

# Female
posterior_mean_S_female <- Cohorts_all_mcmc %>%
  summarize(mean = mean(S_female)) %>%
  pull(mean)
hpd_interval_S_female <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(S_female)), prob = 0.95)[1, ]

# Male
posterior_mean_S_male <- Cohorts_all_mcmc %>%
  summarize(mean = mean(S_male)) %>%
  pull(mean)
hpd_interval_S_male <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(S_male)), prob = 0.95)[1, ]


# Create data frame for values
S_cohorts <- data.frame(
  mean_S = c(posterior_mean_S_female, posterior_mean_S_male),
  hpdi_lower = c(hpd_interval_S_female[1], hpd_interval_S_male[1]),
  hpdi_upper = c(hpd_interval_S_female[2], hpd_interval_S_male[2]),
  sex = c("female", "male"),
  Type = c("cohorts", "cohorts")
)


### Male vs Female: Using all cohorts ###
# posterior prob that S female is > S male

# Count the frequency in the MCMC sample that S_female > S_male = posterior probability that female S is larger than male S
post_prob_M_F <- Cohorts_all_mcmc %>%
  mutate(S_f_minus_m = S_female - S_male) %>%
  summarize(Perc = (sum(S_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_S_Cohorts.txt"))
cat("Posterior probability that S_female > S_male (all cohorts combined):", post_prob_M_F, "\n")
sink()





### Figure S female and male, MA and all cohorts ###

all_data_S <- bind_rows(S_MA, S_cohorts)
  
# Save the combined data frame
save(all_data_S, file = paste0(directory, "S_data_MA_cohorts.RData"))
write.table(all_data_S, file = paste0(directory, "S_data_MA_cohorts.txt"), col.names = T, row.names = F, sep = "\t", quote = F)


S_plot <- ggplot(all_data_S, aes(x = mean_S, y = Type, colour = sex)) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "dashed") +
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
  scale_x_continuous("Selection (S)") +
  annotate("text", x = 1, y = "MA", label = "78%", size = 3) +
  annotate("text", x = 1, y = "cohorts", label = "24%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "right")

S_plot


ggsave(S_plot, file = paste0(directory, "forest_plot_S_MA_cohorts.png"), width = 15, height = 10, unit = "cm")














