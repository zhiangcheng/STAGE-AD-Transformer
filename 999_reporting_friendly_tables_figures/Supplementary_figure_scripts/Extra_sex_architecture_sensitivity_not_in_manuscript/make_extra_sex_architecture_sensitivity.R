# Create publication figure of sensitivity analyses comparing males and females
# SBayesS (heritability, polygenicity, selection paramter) downsampled and across-cohort heterogeneity, h2 varying population prevalence, h2 with varying unscreened controls in males, mixer with downsampling

directory = "/path/10_Publication_figures/Supplementary_figures/Compare_sexes_sensitivity_analyses/"



### Packages ###
library(patchwork)
library(ggforce)      # ggplot circles
library(tidyverse)

# heritability data
load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_MA <- h2_liability_varying_K

load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_u_males.RData")

load("/path/12_SBayesS_h2/GWAS_MA_sumstats_observed_h2.RData")
observed_h2_MA <- observed_h2

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_UKB_full <- h2_liability_varying_K

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_UKB_downsampled <- h2_liability_varying_K


load("/path/12_SBayesS_h2/Cohorts_h2_liability_varying_K.RData")

# polygneicity data

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_polygenicity.RData")
pi_UKB_full <- pi

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_polygenicity.RData")
pi_UKB_downsampled <- pi

load("/path/12_SBayesS_h2/Pi_data_MA_cohorts.RData")
pi_data_cohorts <- all_data_pi %>%
  filter(Type == "cohorts")

# Selection parameter data

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_selection.RData")
S_UKB_full <- S

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_selection.RData")
S_UKB_downsampled <- S

load("/path/12_SBayesS_h2/S_data_MA_cohorts.RData")
S_data_cohorts <- all_data_S %>%
  filter(Type == "cohorts")





### All MCMC samples so can plot distribution
## UK B
# UK B downsampled
UKB_downsample_MCMC_F <- "/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/GWAS_UKB_depression_females_autosomes_and_X_SBayesS.mcmcsamples/CoreParameters.mcmcsamples.txt"

UKB_downsample_MCMC_M <- "/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/GWAS_UKB_depression_males_autosomes_and_X_SBayesS.mcmcsamples/CoreParameters.mcmcsamples.txt"

UKB_downsample_MCMC_F <- read.table(file = UKB_downsample_MCMC_F, header = TRUE, stringsAsFactors = FALSE)

UKB_downsample_MCMC_M <- read.table(file = UKB_downsample_MCMC_M, header = TRUE, stringsAsFactors = FALSE)

UKB_downsample_MCMC_F <- UKB_downsample_MCMC_F %>%
  slice(501:n()) %>%  # remove burnin
  mutate(sex = "female")

UKB_downsample_MCMC_M <- UKB_downsample_MCMC_M %>%
  slice(501:n()) %>%  # remove burnin
  mutate(sex = "male")

UKB_downsample_MCMC_FM <- bind_rows(UKB_downsample_MCMC_F, UKB_downsample_MCMC_M) %>%
  mutate(sample = "UKB_downsampled")

# UK B full sample
UKB_full_MCMC_F <- "/path/12_SBayesS_h2/UKB_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt"

UKB_full_MCMC_M <- "/path/12_SBayesS_h2/UKB_AD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt"

UKB_full_MCMC_F <- read.table(file = UKB_full_MCMC_F, header = TRUE, stringsAsFactors = FALSE)

UKB_full_MCMC_M <- read.table(file = UKB_full_MCMC_M, header = TRUE, stringsAsFactors = FALSE)

UKB_full_MCMC_F <- UKB_full_MCMC_F %>%
  slice(501:n()) %>%  # remove burnin
  mutate(sex = "female")

UKB_full_MCMC_M <- UKB_full_MCMC_M %>%
  slice(501:n()) %>%  # remove burnin
  mutate(sex = "male")

UKB_full_MCMC_FM <- bind_rows(UKB_full_MCMC_F, UKB_full_MCMC_M) %>%
  mutate(sample = "UKB_full")


## Across cohort heterogeneity: All MCMC samples from all cohorts combined
SBayesS_directory = "/path/12_SBayesS_h2/"

AGDS_females_mcmc <- paste0(SBayesS_directory, "AGDS_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
AGDS_males_mcmc <- paste0(SBayesS_directory, "AGDS_AD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

AllOfUs_females_mcmc <- paste0(SBayesS_directory, "AllOfUs_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
AllOfUs_males_mcmc <- paste0(SBayesS_directory, "AllOfUs_AD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

Bionic_females_mcmc <- paste0(SBayesS_directory, "Bionic_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
Bionic_males_mcmc <- paste0(SBayesS_directory, "Bionic_AD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

Blokland21_females_mcmc <- paste0(SBayesS_directory, "Blokland21_AD_female_sumstats_SBayesS_DENTISTremoved_redo.mcmcsamples/CoreParameters.mcmcsamples.txt")
Blokland21_males_mcmc <- paste0(SBayesS_directory, "Blokland21_AD_male_sumstats_trimmedN_30perc_SBayesS.mcmcsamples/CoreParameters.mcmcsamples.txt")

GLAD_females_mcmc <- paste0(SBayesS_directory, "GLAD_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
GLAD_males_mcmc <- paste0(SBayesS_directory, "GLAD_AD_male_sumstats_SBayesS_h006_pi001.mcmcsamples/CoreParameters.mcmcsamples.txt")

UKB_females_mcmc <- paste0(SBayesS_directory, "UKB_AD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
UKB_males_mcmc <- paste0(SBayesS_directory, "UKB_AD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

read_and_rename_mcmc <- function(file_path, Sex, Study) {
  read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    slice(501:n()) %>%
    mutate(sex = Sex,
           study = Study)
  
}

# Process each cohort
AGDS_F_mcmc <- read_and_rename_mcmc(AGDS_females_mcmc, "female", "AGDS")
AGDS_M_mcmc <- read_and_rename_mcmc(AGDS_females_mcmc, "male", "AGDS")

AllOfUs_F_mcmc <- read_and_rename_mcmc(AllOfUs_females_mcmc, "female", "AllOfUs")
AllOfUs_M_mcmc <- read_and_rename_mcmc(AllOfUs_females_mcmc, "male", "AllOfUs")

Bionic_F_mcmc <- read_and_rename_mcmc(Bionic_females_mcmc, "female", "Bionic")
Bionic_M_mcmc <- read_and_rename_mcmc(Bionic_females_mcmc, "male", "Bionic")

Blokland21_F_mcmc <- read_and_rename_mcmc(Blokland21_females_mcmc, "female", "Blokland")
Blokland21_M_mcmc <- read_and_rename_mcmc(Blokland21_females_mcmc, "male", "Blokland")

GLAD_F_mcmc <- read_and_rename_mcmc(GLAD_females_mcmc, "female", "GLAD")
GLAD_M_mcmc <- read_and_rename_mcmc(GLAD_females_mcmc, "male", "GLAD")

UKB_F_mcmc <- read_and_rename_mcmc(UKB_females_mcmc, "female", "UKB")
UKB_M_mcmc <- read_and_rename_mcmc(UKB_females_mcmc, "male", "UKB")

## Sex-stratified meta-analysis
SBayes_directory = "/path/12_SBayesS_h2/"
MA_females_mcmc <- paste0(SBayes_directory, "Metaanalysis_AD_female_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
MA_males_mcmc <- paste0(SBayes_directory, "Metaanalysis_AD_male_AllCohorts_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")

MA_females_mcmc_df <- read.table(MA_females_mcmc, header = TRUE, stringsAsFactors = FALSE)
MA_males_mcmc_df <- read.table(MA_males_mcmc, header = TRUE, stringsAsFactors = FALSE)

MA_females_mcmc_df <- MA_females_mcmc_df %>%
  mutate(sex = "female") %>% 
  slice(501:n())

MA_males_mcmc_df <- MA_males_mcmc_df %>%
  mutate(sex = "male") %>% 
  slice(501:n())

MA_all_mcmc <- bind_rows(MA_females_mcmc_df, MA_males_mcmc_df) 


### Calculate h2 on liability scale ###

# Function to convert h2 to liability scale
# Lee et al 2011 AJHG 's method to convert the heritability estimate and standard error at the observed scale to those at the liability scale obs is the estimate (or SE) at the observed scale K is the population prevalence P is the sample prevalence

mapToLiabilityScale = function(obs, K, P){
  z = dnorm(qnorm(1-K))
  lia = obs * (K*(1-K)/z^2) * (K*(1-K)/(P*(1-P)))
  return(lia)
}

# Convert all hsq values in MCMC samples to liability scale

## UK B
UKB_full_MCMC_FM <- UKB_full_MCMC_FM %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.20, P = 0.4647),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.10, P = 0.2857))
  )

UKB_downsample_MCMC_FM <- UKB_downsample_MCMC_FM %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.20, P = 0.298),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.10, P = 0.298))
  )

UKB_MCMC <- bind_rows(UKB_downsample_MCMC_FM, UKB_full_MCMC_FM)

## Each cohort individually (to look at across cohort heterogeneity)
calculate_liability <- function(df, cohort, P_values) {
  df %>%
    mutate(
      hsq_l = case_when(
        sex == "female" ~ mapToLiabilityScale(hsq, K = 0.20, P = P_values$female),
        sex == "male" ~ mapToLiabilityScale(hsq, K = 0.10, P = P_values$male)
    ))
}

P_values <- list(
  AGDS = list(female = 0.5928, male = 0.3247),
  AllOfUs = list(female = 0.3645, male = 0.2371),
  Bionic = list(female = 0.2841, male = 0.1813),
  Blokland21 = list(female = 0.4703, male = 0.6484),
  GLAD = list(female = 0.8312, male = 0.6062),
  UKB = list(female = 0.4647, male = 0.2857)
)

# Apply the function to each cohort's data
AGDS_F_mcmc <- calculate_liability(AGDS_F_mcmc, "AGDS", P_values$AGDS)
AGDS_M_mcmc <- calculate_liability(AGDS_M_mcmc, "AGDS", P_values$AGDS)
AllOfUs_F_mcmc <- calculate_liability(AllOfUs_F_mcmc, "AllOfUs", P_values$AllOfUs)
AllOfUs_M_mcmc <- calculate_liability(AllOfUs_M_mcmc, "AllOfUs", P_values$AllOfUs)
Bionic_F_mcmc <- calculate_liability(Bionic_F_mcmc, "Bionic", P_values$Bionic)
Bionic_M_mcmc <- calculate_liability(Bionic_M_mcmc, "Bionic", P_values$Bionic)
Blokland21_F_mcmc <- calculate_liability(Blokland21_F_mcmc, "Blokland21", P_values$Blokland21)
Blokland21_M_mcmc <- calculate_liability(Blokland21_M_mcmc, "Blokland21", P_values$Blokland21)
GLAD_F_mcmc <- calculate_liability(GLAD_F_mcmc, "GLAD", P_values$GLAD)
GLAD_M_mcmc <- calculate_liability(GLAD_M_mcmc, "GLAD", P_values$GLAD)
UKB_F_mcmc <- calculate_liability(UKB_F_mcmc, "UKB", P_values$UKB)
UKB_M_mcmc <- calculate_liability(UKB_M_mcmc, "UKB", P_values$UKB)

# Combine into one df
Cohorts_MCMC <- bind_rows(
  AGDS_F_mcmc,
  AGDS_M_mcmc,
  AllOfUs_F_mcmc,
  AllOfUs_M_mcmc,
  Bionic_F_mcmc,
  Bionic_M_mcmc,
  Blokland21_F_mcmc,
  Blokland21_M_mcmc,
  GLAD_F_mcmc,
  GLAD_M_mcmc,
  UKB_F_mcmc,
  UKB_M_mcmc
)

## Sex-stratified meta-analysis GWAS, use varying K to convert to h2 liability
MA_all_mcmc_005 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.05, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.05, P = 0.3290),
  ),
  K = 0.05)

MA_all_mcmc_010 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.10, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.10, P = 0.3290),
  ),
  K = 0.10)

MA_all_mcmc_015 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.15, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.15, P = 0.3290),
  ),
  K = 0.15)

MA_all_mcmc_020 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.20, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.20, P = 0.3290),
  ),
  K = 0.20)

MA_all_mcmc_025 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScale(hsq, K = 0.25, P = 0.4499),
    sex == "male" ~ mapToLiabilityScale(hsq, K = 0.25, P = 0.3290),
  ),
  K = 0.25)

MA_MCMC <- bind_rows(
  MA_all_mcmc_005,
  MA_all_mcmc_010,
  MA_all_mcmc_015,
  MA_all_mcmc_020,
  MA_all_mcmc_025
) %>% 
  mutate(sex = factor(sex, levels = c("male", "female")))

## Sex-stratified meta-analysis GWAS, use varying K and u for males to convert to h2 liability

# Create function to convert to liability scale taking into account unscreened controls - Peyrot et al., 2016
# F = proportion of falsely classified control subjects (F = N(false controls) / N(total controls))
# u = proportion of unscreened control subjects (F = K*u)

mapToLiabilityScaleUnscreened = function(obs, K, P, u){
  z = dnorm(qnorm(1-K))
  lia = obs * (K^2*(1-K)^2)/(P*(1-P)*(1-K*u)^2*z^2)
  return(lia)
}

MA_all_mcmc_F <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "female" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.2, P = 0.4499, u = 0)
  ),
  K = 0.20,
  u = 0.0)

MA_all_mcmc_1 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.2, P = 0.3290, u = 1)
  ),
  K = 0.20,
  u = 1.0)

MA_all_mcmc_2 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.19, P = 0.3290, u = 0.9)
  ),
  K = 0.19,
  u = 0.9)

MA_all_mcmc_3 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.18, P = 0.3290, u = 0.8)
  ),
  K = 0.18,
  u = 0.8)

MA_all_mcmc_4 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.17, P = 0.3290, u = 0.7)
  ),
  K = 0.17,
  u = 0.7)

MA_all_mcmc_5 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.16, P = 0.3290, u = 0.6)
  ),
  K = 0.16,
  u = 0.6)

MA_all_mcmc_6 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.15, P = 0.3290, u = 0.5)
  ),
  K = 0.15,
  u = 0.5)

MA_all_mcmc_7 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.14, P = 0.3290, u = 0.4)
  ),
  K = 0.14,
  u = 0.4)

MA_all_mcmc_8 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.13, P = 0.3290, u = 0.3)
  ),
  K = 0.13,
  u = 0.3)

MA_all_mcmc_9 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.12, P = 0.3290, u = 0.2)
  ),
  K = 0.12,
  u = 0.2)

MA_all_mcmc_10 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.11, P = 0.3290, u = 0.1)
  ),
  K = 0.11,
  u = 0.1)

MA_all_mcmc_11 <- MA_all_mcmc %>%
  mutate(hsq_l = case_when(
    sex == "male" ~ mapToLiabilityScaleUnscreened(hsq, K = 0.10, P = 0.3290, u = 0)
  ),
  K = 0.10,
  u = 0.0)

MA_MCMC_varyingu <- bind_rows(
  MA_all_mcmc_F,
  MA_all_mcmc_1,
  MA_all_mcmc_2,
  MA_all_mcmc_3,
  MA_all_mcmc_4,
  MA_all_mcmc_5,
  MA_all_mcmc_6,
  MA_all_mcmc_7,
  MA_all_mcmc_8,
  MA_all_mcmc_9,
  MA_all_mcmc_10,
  MA_all_mcmc_11
) %>%
  mutate(sex = factor(sex, levels = c("male", "female")))

### SNP-based heritability Plot: Liability scale with UKB full and downsampled ###

h2_liability_UKB_full <- h2_liability_varying_K_UKB_full %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(cohort = "UKB_full") %>%
  select(-K)

h2_liability_UKB_downsampled <- h2_liability_varying_K_UKB_downsampled %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(cohort = "UKB_downsampled") %>%
  select(-K)

h2_liability_UKB <- h2_liability_UKB_full %>%
  bind_rows(h2_liability_UKB_downsampled)

h2_plot_sexes_UKB <- ggplot(h2_liability_UKB, aes(x = cohort, y = mean_hsq * 100, colour = sex)) +
  geom_violin(data = UKB_MCMC, aes(x = sample, y = hsq_l * 100, colour = sex, fill = sex),
              position = position_dodge(width = -0.6),
              inherit.aes = F,
              show.legend = F) +
  geom_point(position = position_dodge(width = -0.6)) +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100),
                 position = position_dodge(width = -0.6)) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(9, 23),
                     breaks = seq(10, 22, 2)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  annotate("text", x = 1, y = 22, label = "P(F > M) \n= 99%", size = 3) +
  annotate("text", x = 2, y = 22, label = "P(F > M) \n= 100%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
h2_plot_sexes_UKB


### Polygenicity Plot: UKB full and downsampled ###

pi_UKB <- (pi_UKB_full %>% mutate(cohort = "UKB_full")) %>%
  bind_rows(pi_UKB_downsampled %>% mutate(cohort = "UKB_downsampled"))


pi_plot_sexes_UKB <- ggplot(pi_UKB, aes(x = cohort, y = mean_pi, colour = sex)) +
  geom_violin(data = UKB_MCMC, aes(x = sample, y = Pi, colour = sex, fill = sex),
              position = position_dodge(width = -0.6),
              inherit.aes = F,
              show.legend = F) +
  geom_point(position = position_dodge(width = -0.6)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.6)) +
  scale_y_continuous("Polygenicity (\u03c0)",
                     limits = c(0, 0.028),
                     breaks = seq(0, 0.024, 0.004)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  annotate("text", x = 1, y = 0.026, label = "P(F > M) \n= 87%", size = 3) +
  annotate("text", x = 2, y = 0.026, label = "P(F > M) \n= 29%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
pi_plot_sexes_UKB


### Selection Parameter Plot: UKB full and downsampled ###

S_UKB <- (S_UKB_full %>% mutate(cohort = "UKB_full")) %>%
  bind_rows(S_UKB_downsampled %>% mutate(cohort = "UKB_downsampled"))


S_plot_sexes_UKB <- ggplot(S_UKB, aes(x = cohort, y = mean_S, colour = sex)) +
  geom_violin(data = UKB_MCMC, aes(x = sample, y = S, colour = sex, fill = sex),
              position = position_dodge(width = -0.6),
              inherit.aes = F,
              show.legend = F) +
  geom_point(position = position_dodge(width = -0.6)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.6)) +
  scale_y_continuous(expression("Selection parameter (" * italic(S) * ")"),
                     limits = c(-0.9, 0.4),
                     breaks = seq(-0.8, 0.2, 0.2)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  annotate("text", x = 1, y = 0.3, label = "P(F > M) \n= 79%", size = 3) +
  annotate("text", x = 2, y = 0.3, label = "P(F > M) \n= 25%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
S_plot_sexes_UKB




### SNP-based heritability Plot: Across cohort heterogeneity ###

h2_liability <- h2_liability_varying_K_cohorts %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(scale = "Liability") %>%
  select(-K)


h2_plot_sexes_cohorts <- ggplot(h2_liability, aes(x = sex, y = mean_hsq * 100, colour = sex)) +
  geom_violin(data = Cohorts_MCMC, aes(x = sex, y = hsq * 100, colour = sex, fill = sex),
             inherit.aes = F,
             show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100)) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(0, 51)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_manual(
    values = c('female' = '#FDE725FF',
               'male' = '#440154FF'),
    breaks = c("male", "female"),
    labels = c('male' = "Male",
               'female' = "Female")
  ) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  annotate("text", x = 1.5, y = 51, label = "P(F > M) = 93%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_plot_sexes_cohorts


### Polygenicity Plot: Across-cohort heterogeneity ###


pi_plot_sexes_cohorts <- ggplot(pi_data_cohorts, aes(x = sex, y = mean_pi, colour = sex)) +
  geom_violin(data = Cohorts_MCMC, aes(x = sex, y = Pi, colour = sex, fill = sex),
              inherit.aes = F,
              show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_y_continuous("Polygenicity (\u03c0)",
                     limits = c(0,0.045)) +
  annotate("text", x = 1.5, y = 0.044, label = "P(F > M) = 79%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

pi_plot_sexes_cohorts

### Selection Plot: Across-cohort heterogeneity ###

S_plot_sexes_cohorts <- ggplot(S_data_cohorts, aes(x = sex, y = mean_S, colour = sex)) +
  geom_violin(data = Cohorts_MCMC, aes(x = sex, y = S, colour = sex, fill = sex),
              inherit.aes = F,
              show.legend = F) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  scale_fill_viridis_d(end = 0, begin = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_y_continuous(expression("Selection parameter (" * italic(S) * ")"),
                     limits = c(-1.6, 2.5),
                     breaks = seq(-1.5, 2.5, 0.5)) +
  annotate("text", x = 1.5, y = 2.5, label = "P(F > M) = 24%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

S_plot_sexes_cohorts


### SNP-based heritability Plot: Multiple K (SBayesS) ###

# Using results from SBayesS on GWAS sex stratified meta-analysis results

h2_liability_varying_K_MA_plot <- h2_liability_varying_K_MA %>%
  mutate(sex = factor(sex, levels = c("male", "female")))

h2_varying_K_plot <- ggplot(h2_liability_varying_K_MA_plot, aes(x = as.factor(K), y = mean_hsq*100, colour = sex)) +
  geom_vline(xintercept = 2,
             colour = "grey30",
             linetype = 'dashed') +
  geom_vline(xintercept = 4,
             colour = "grey30",
             linetype = 'dashed') +
  geom_violin(data = MA_MCMC, aes(x = as.factor(K), y = hsq_l*100, colour = as.factor(sex), fill = as.factor(sex)),
              position = position_dodge(width = 0.6),
              inherit.aes = F,
              show.legend = F) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.6)) +
  scale_colour_viridis_d(begin = 0, end = 1) +
  scale_fill_viridis_d(begin = 0, end = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_x_discrete("Population Prevalence") +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(6, 14),
                     breaks = seq(6, 14, 2)) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_K_plot


### SNP-based heritability Plot: Varying u (SBayesS) ###

h2_liability_varying_u <- h2_liability_varying_u_males %>%
  mutate(sex = "male") %>%
  slice_tail(n = 11) %>%
  bind_rows(h2_liability_varying_K_MA %>% filter(sex == "female" & K == 0.20) %>% mutate(u = 0)) %>%
  mutate(sex = factor(sex, levels = c("male", "female")))

h2_varying_u_plot <- ggplot(h2_liability_varying_u, aes(x = as.factor(u), y = mean_hsq*100, colour = sex)) +
  geom_violin(data = MA_MCMC_varyingu, aes(x = as.factor(u), y = hsq_l*100, colour = sex, fill = sex),
              position = position_dodge(width = 0.4),
              inherit.aes = F,
              show.legend = F) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.4)) +
  scale_colour_viridis_d(begin = 0, end = 1) +
  scale_fill_viridis_d(begin = 0, end = 1,
                       alpha = 0.1) +
  guides(color = "none") +
  scale_x_discrete("Proportion of unscreened control subjects",
                   breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expression("Liability " * italic(h)^2 * " (%)"),
                     limits = c(8, 20),
                     breaks = seq(8, 20, 2)) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_u_plot



### Mixer Venn: Full UKB ###

# Data: sizes of each group
one <- 3725  # female exclusive
two <- 104   # male exclusive
overlap <- 7516  # Shared between female and male

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
rm(pi)
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_UKB_full <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  guides(fill = "none") +
  coord_fixed(expand = TRUE) +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = mean(circle_data$x) + 0.6,
           y = circle_data$y[1],
           label = overlap,
           size = 3, color = "black", family = "Calibri") +
  ggtitle("UK B Full") +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


mixer_venn_UKB_full



### Mixer Venn: Downsampled UKB ###

# Data: sizes of each group
one <- 5769  # female exclusive
two <- 132   # male exclusive
overlap <- 6645  # Shared between female and male

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
rm(pi)
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("AD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_UKB_down <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  guides(fill = "none") +
  coord_fixed(expand = TRUE) +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = mean(circle_data$x) + 0.6,
           y = circle_data$y[1],
           label = overlap,
           size = 3, color = "black", family = "Calibri") +
  ggtitle("UK B Downsampled") +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


mixer_venn_UKB_down


### Create Panelled Figure ###

figure <- (guide_area() /
             ((h2_plot_sexes_UKB | pi_plot_sexes_UKB | S_plot_sexes_UKB ) + plot_layout(widths = c(1, 1, 1))) /
             ((h2_plot_sexes_cohorts | pi_plot_sexes_cohorts | S_plot_sexes_cohorts) + plot_layout(widths = c(1, 1, 1))) /
             ((h2_varying_K_plot | h2_varying_u_plot) + plot_layout(widths = c(1, 1.5))) /
             ((mixer_venn_UKB_full | mixer_venn_UKB_down))) +
  plot_layout(heights = c(0.2, 1, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "H", "I", ""))) &
  theme(legend.position = "top")

figure


outfile <- paste(directory, "Compare_sexes_sensitivity_analyses.png", sep="")
ggsave(figure, width = 21, height = 29, unit = "cm", file = outfile)
