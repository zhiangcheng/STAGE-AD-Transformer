# Compare male and female h2 from meta-analysed sumstats and after combining all per cohort h2 estimates
# Using lifetime pop prev of 20% for females and 10% for males

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
Blokland21_males_mcmc <- paste0(directory, "Blokland21_MDD_male_sumstats_trimmedN_30perc_SBayesS.mcmcsamples/CoreParameters.mcmcsamples.txt")

GLAD_females_mcmc <- paste0(directory, "GLAD_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
GLAD_males_mcmc <- paste0(directory, "GLAD_MDD_male_sumstats_SBayesS_h006_pi001.mcmcsamples/CoreParameters.mcmcsamples.txt")

UKB_females_mcmc <- paste0(directory, "UKB_MDD_female_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")
UKB_males_mcmc <- paste0(directory, "UKB_MDD_male_sumstats_SBayesS_redomcmc.mcmcsamples/CoreParameters.mcmcsamples.txt")



### Packages ###

library(metafor)      # meta-analysis
library(coda)         # HPDI
library(tidyverse)



### Load Data ###

# Load mcmc reports for sex-stratified MA sumstats

MA_females_mcmc_df <- read.table(MA_females_mcmc, header = TRUE, stringsAsFactors = FALSE)
MA_males_mcmc_df <- read.table(MA_males_mcmc, header = TRUE, stringsAsFactors = FALSE)

MA_females_mcmc_df <- MA_females_mcmc_df %>%
  rename_with(~ paste0(., "_female"))

MA_males_mcmc_df <- MA_males_mcmc_df %>%
  rename_with(~ paste0(., "_male"))

MA_all_mcmc <- bind_cols(MA_females_mcmc_df, MA_males_mcmc_df)



### Function to convert h2 to liability scale ###

# Lee et al 2011 AJHG 's method to convert the heritability estimate and standard error at the observed scale to those at the liability scale obs is the estimate (or SE) at the observed scale K is the population prevalence P is the sample prevalence

mapToLiabilityScale = function(obs, K, P){
  z = dnorm(qnorm(1-K))
  lia = obs * (K*(1-K)/z^2) * (K*(1-K)/(P*(1-P)))
  return(lia)
}
# mapToLiabilityScaleVar = function(obs, K, P){
#   z = dnorm(qnorm(1-K))
#   lia = obs * ((K*(1-K)/z^2) * (K*(1-K)/(P*(1-P))))^2
#   return(lia)
# }


# Create function to convert to liability scale taking into account unscreened controls - Peyrot et al., 2016
# F = proportion of falsely classified control subjects (F = N(false controls) / N(total controls))
# u = proportion of unscreened control subjects (F = K*u)

mapToLiabilityScaleUnscreened = function(obs, K, P, u){
  z = dnorm(qnorm(1-K))
  lia = obs * (K^2*(1-K)^2)/(P*(1-P)*(1-K*u)^2*z^2)
  return(lia)
}


### Convert all hsq values in MCMC samples to liability scale (sex-stratified MA sumstats) ###
MA_all_mcmc <- MA_all_mcmc %>%
  mutate(hsq_F_L25 = mapToLiabilityScale(hsq_female, K = 0.25, P = 0.4499),
         hsq_F_L20 = mapToLiabilityScale(hsq_female, K = 0.20, P = 0.4499),
         hsq_F_L15 = mapToLiabilityScale(hsq_female, K = 0.15, P = 0.4499),
         hsq_F_L10 = mapToLiabilityScale(hsq_female, K = 0.10, P = 0.4499),
         hsq_F_L5 = mapToLiabilityScale(hsq_female, K = 0.05, P = 0.4499),
         
         hsq_M_L25 = mapToLiabilityScale(hsq_male, K = 0.25, P = 0.3290),
         hsq_M_L20 = mapToLiabilityScale(hsq_male, K = 0.20, P = 0.3290),
         hsq_M_L15 = mapToLiabilityScale(hsq_male, K = 0.15, P = 0.3290),
         hsq_M_L10 = mapToLiabilityScale(hsq_male, K = 0.10, P = 0.3290),
         hsq_M_L5 = mapToLiabilityScale(hsq_male, K = 0.05, P = 0.3290),
         
         hsq_M_L10_unscr10 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 1),
         hsq_M_L10_unscr9 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.9),
         hsq_M_L10_unscr8 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.8),
         hsq_M_L10_unscr7 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.7),
         hsq_M_L10_unscr6 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.6),
         hsq_M_L10_unscr5 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.5),
         hsq_M_L10_unscr4 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.4),
         hsq_M_L10_unscr3 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.3),
         hsq_M_L10_unscr2 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.2),
         hsq_M_L10_unscr1 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0.1),
         hsq_M_L10_unscr0 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0),
         
         hsq_M_L20_unscr10 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.20, P = 0.3290, u = 1),
         hsq_M_L19_unscr9 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.19, P = 0.3290, u = 0.9),
         hsq_M_L18_unscr8 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.18, P = 0.3290, u = 0.8),
         hsq_M_L17_unscr7 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.17, P = 0.3290, u = 0.7),
         hsq_M_L16_unscr6 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.16, P = 0.3290, u = 0.6),
         hsq_M_L15_unscr5 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.15, P = 0.3290, u = 0.5),
         hsq_M_L14_unscr4 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.14, P = 0.3290, u = 0.4),
         hsq_M_L13_unscr3 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.13, P = 0.3290, u = 0.3),
         hsq_M_L12_unscr2 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.12, P = 0.3290, u = 0.2),
         hsq_M_L11_unscr1 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.11, P = 0.3290, u = 0.1),
         hsq_M_L10_unscr0 = mapToLiabilityScaleUnscreened(hsq_male, K = 0.10, P = 0.3290, u = 0))




### Liability h2 for a range of K ###
# Create a list of column names, sex labels, and K values for each column
liability_columns <- c("hsq_F_L25", "hsq_F_L20", "hsq_F_L15", "hsq_F_L10", "hsq_F_L5",
                       "hsq_M_L25", "hsq_M_L20", "hsq_M_L15", "hsq_M_L10", "hsq_M_L5")

sex_labels <- c(rep("female", 5), rep("male", 5))
K_values <- c(0.25, 0.20, 0.15, 0.10, 0.05, 0.25, 0.20, 0.15, 0.10, 0.05)


# Initialize an empty list to store results for liability scale
results_list <- list()

# Iterate over each column and calculate posterior mean and 95% HPD interval for liability scale
for (i in seq_along(liability_columns)) {
  column_name <- liability_columns[i]
  sex <- sex_labels[i]
  K <- K_values[i]
  
  # Calculate the posterior mean for liability scale
  # 2,500 data points = chain length of 25,000 with thinning of 10
  # burnin = 5,000 (so 500 data point when thinning of 10)
  posterior_mean_liability <- MA_all_mcmc %>%
    slice(501:n()) %>%
    summarize(mean = mean(.data[[column_name]])) %>%
    pull(mean)
  
  # Convert to MCMC object and calculate 95% HPD interval for liability scale
  mcmc_obj_liability <- as.mcmc(MA_all_mcmc %>% slice(501:n()) %>% pull(.data[[column_name]]))
  hpd_interval_liability <- HPDinterval(mcmc_obj_liability, prob = 0.95)[1, ]
  
  # Store liability scale results in the list
  results_list[[i]] <- data.frame(
    mean_hsq = posterior_mean_liability,
    hpdi_lower = hpd_interval_liability[1],
    hpdi_upper = hpd_interval_liability[2],
    sex = sex,
    K = K
  )
}

# Combine all liability-scale results into a single data frame
h2_liability_varying_K <- do.call(rbind, results_list)



### Observed-scale heritability for female and male ###
# Female observed scale
observed_female <- MA_all_mcmc %>%
  slice(501:n()) %>%
  summarize(mean = mean(hsq_female)) %>%
  pull(mean)

observed_female_hpd <- HPDinterval(as.mcmc(MA_all_mcmc %>% slice(501:n()) %>% pull(hsq_female)), prob = 0.95)[1, ]

# Male observed scale
observed_male <- MA_all_mcmc %>%
  slice(501:n()) %>%
  summarize(mean = mean(hsq_male)) %>%
  pull(mean)

observed_male_hpd <- HPDinterval(as.mcmc(MA_all_mcmc %>% slice(501:n()) %>% pull(hsq_male)), prob = 0.95)[1, ]

# Create data frame for observed values and bind to liability scale results
observed_h2 <- data.frame(
  mean_hsq = c(observed_female, observed_male),
  hpdi_lower = c(observed_female_hpd[1], observed_male_hpd[1]),
  hpdi_upper = c(observed_female_hpd[2], observed_male_hpd[2]),
  sex = c("female", "male"),
  K = "observed"
)


### Liability scale heritability for males using a range of u and 10% pop prev or changing pop prev with u ###
# Create a list of column names, sex labels, and K values for each column

liability_columns <- c("hsq_M_L10_unscr10", "hsq_M_L10_unscr9", "hsq_M_L10_unscr8", "hsq_M_L10_unscr7", "hsq_M_L10_unscr6", "hsq_M_L10_unscr5", "hsq_M_L10_unscr4", "hsq_M_L10_unscr3", "hsq_M_L10_unscr2", "hsq_M_L10_unscr1", "hsq_M_L10_unscr0", "hsq_M_L20_unscr10", "hsq_M_L19_unscr9", "hsq_M_L18_unscr8", "hsq_M_L17_unscr7", "hsq_M_L16_unscr6", "hsq_M_L15_unscr5", "hsq_M_L14_unscr4", "hsq_M_L13_unscr3", "hsq_M_L12_unscr2", "hsq_M_L11_unscr1", "hsq_M_L10_unscr0")

u_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)
K_values <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1)

# Initialize an empty list to store results for liability scale
results_list <- list()


# Iterate over each column and calculate posterior mean and 95% HPD interval for liability scale
for (i in seq_along(liability_columns)) {
  column_name <- liability_columns[i]
  u <- u_values[i]
  K <- K_values[i]
  
  # Calculate the posterior mean for liability scale
  # 2,500 data points = chain length of 25,000 with thinning of 10
  # burnin = 5,000 (so 500 data point when thinning of 10)
  posterior_mean_liability <- MA_all_mcmc %>%
    slice(501:n()) %>%
    summarize(mean = mean(.data[[column_name]])) %>%
    pull(mean)
  
  # Convert to MCMC object and calculate 95% HPD interval for liability scale
  mcmc_obj_liability <- as.mcmc(MA_all_mcmc %>% slice(501:n()) %>% pull(.data[[column_name]]))
  hpd_interval_liability <- HPDinterval(mcmc_obj_liability, prob = 0.95)[1, ]
  
  # Store liability scale results in the list
  results_list[[i]] <- data.frame(
    mean_hsq = posterior_mean_liability,
    hpdi_lower = hpd_interval_liability[1],
    hpdi_upper = hpd_interval_liability[2],
    K = K,
    u = u
  )
}

# Combine all liability-scale results into a single data frame
h2_liability_varying_u_males <- do.call(rbind, results_list)




# Save the combined data frame
save(h2_liability_varying_K, file = paste0(directory, "GWAS_MA_sumstats_h2_liability_varying_K.RData"))
save(observed_h2, file = paste0(directory, "GWAS_MA_sumstats_observed_h2.RData"))
save(h2_liability_varying_u_males, file = paste0(directory, "GWAS_MA_sumstats_h2_liability_varying_u_males.RData"))

write.table(h2_liability_varying_K, file = paste0(directory, "GWAS_MA_sumstats_h2_liability_varying_K.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(observed_h2, file = paste0(directory, "GWAS_MA_sumstats_observed_h2.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(h2_liability_varying_u_males, file = paste0(directory, "GWAS_MA_sumstats_h2_liability_varying_u_males.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



### Forest Plot liability h2 with varying K ###

h2_varying_K_plot <- ggplot((h2_liability_varying_K %>% filter(K != 0.0563 & K != 0.1018)), aes(x = K, y = mean_hsq*100, colour = sex)) +
  geom_vline(xintercept = 0.10,
             colour = "grey30",
             linetype = 'dashed') +
  geom_vline(xintercept = 0.20,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 1, end = 0) +
  scale_x_continuous("Population Prevalence",
                     breaks = c(0.05, 0.1, 0.15, 0.20, 0.25),
                     labels = c(0.05, 0.1, 0.15, 0.20, 0.25)) +
  scale_y_continuous(expression(italic(h)^2 * " (%)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_K_plot


ggsave(h2_varying_K_plot, file = paste0(directory, "h2_liability_M_F_varying_K_plot_SBayesS_GWAS_MA_sumstats.png", sep = ""), width = 10, height = 10, unit = "cm")



### Forest Plot liability h2 for males with varying F (K = 10% for all) ###

h2_liability_varying_u <- h2_liability_varying_u_males %>%
  mutate(sex = "male") %>%
  filter(K == 0.1) %>% 
  bind_rows(h2_liability_varying_K %>% filter(sex == "female" & K == 0.20) %>% mutate(u = 0))

h2_varying_u_plot <- ggplot(h2_liability_varying_u, aes(x = u, y = mean_hsq*100, colour = sex)) +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 1, end = 0) +
  scale_x_continuous("Proportion of unscreened control subjects",
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(expression("Liability " * italic(h)^2 * " (%)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_u_plot


ggsave(h2_varying_u_plot, file = paste0(directory, "h2_liability_M_varying_u_plot_K10_SBayesS_GWAS_MA_sumstats.png", sep = ""), width = 10, height = 10, unit = "cm")


### Forest Plot liability h2 for males with varying F and varying K ###

h2_liability_varying_u <- h2_liability_varying_u_males %>%
  mutate(sex = "male") %>%
  slice_tail(n = 11) %>% 
  bind_rows(h2_liability_varying_K %>% filter(sex == "female" & K == 0.20) %>% mutate(u = 0))

h2_varying_u_plot <- ggplot(h2_liability_varying_u, aes(x = u, y = mean_hsq*100, colour = sex)) +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 1, end = 0) +
  scale_x_continuous("Proportion of unscreened control subjects",
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(expression("Liability " * italic(h)^2 * " (%)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_u_plot


ggsave(h2_varying_u_plot, file = paste0(directory, "h2_liability_M_varying_u_varying_K_plot_SBayesS_GWAS_MA_sumstats.png", sep = ""), width = 10, height = 10, unit = "cm")




### Male vs Female: GWAS M-A sumstats ###
# posterior prob that h2 female is > h2 male

# Count the frequency in the MCMC sample that h2_female > h2_male = posterior probability that female h2 is larger than male h2
post_prob_M_F <- MA_all_mcmc %>%
  slice(501:n()) %>%
  mutate(hsq_f_minus_m = hsq_female - hsq_male) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

post_prob_M_F_liability <- MA_all_mcmc %>%
  slice(501:n()) %>%
  mutate(hsq_f_minus_m = hsq_F_L20 - hsq_M_L10) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)



liability_columns <- c("hsq_M_L10_unscr10", "hsq_M_L10_unscr9", "hsq_M_L10_unscr8", "hsq_M_L10_unscr7", "hsq_M_L10_unscr6", "hsq_M_L10_unscr5", "hsq_M_L10_unscr4", "hsq_M_L10_unscr3", "hsq_M_L10_unscr2", "hsq_M_L10_unscr1", "hsq_M_L10_unscr0", "hsq_M_L20_unscr10", "hsq_M_L19_unscr9", "hsq_M_L18_unscr8", "hsq_M_L17_unscr7", "hsq_M_L16_unscr6", "hsq_M_L15_unscr5", "hsq_M_L14_unscr4", "hsq_M_L13_unscr3", "hsq_M_L12_unscr2", "hsq_M_L11_unscr1", "hsq_M_L10_unscr0")

results_list <- list()

for (i in seq_along(liability_columns)) {
  column_name <- liability_columns[i]
  
  post_prob_M_F <- MA_all_mcmc %>%
    slice(501:n()) %>%
    mutate(hsq_f_minus_m = hsq_F_L20 - .data[[column_name]]) %>%
    summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
    mutate(Perc = round(Perc, 2)) %>%
    pull(Perc)
  
  results_list[[i]] <- data.frame(
    male_liability_h2 = column_name,
    post_prob_F_larger_M = post_prob_M_F
  )
}

PP_h2_F_M_varying_u_males <- do.call(rbind, results_list)

write.table(PP_h2_F_M_varying_u_males, file = paste0(directory, "male_vs_female_h2_GWAS_MA_sumstats_varying_u.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Print the result
sink(paste0(directory, "/male_vs_female_h2_GWAS_MA_sumstats.txt"))
cat("Posterior probability that h2_female (observed scale) > h2_male (observed scale):", post_prob_M_F, "\n")
cat("Posterior probability that h2_female (liability scale with K = 0.20) > h2_male (liability scale with K = 0.10):", post_prob_M_F_liability, "\n")
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





### Convert all hsq values in MCMC samples to liability scale ###

calculate_liability <- function(df, cohort, P_values) {
  df %>%
    mutate(
      hsq_F_L25 = mapToLiabilityScale(hsq_female, K = 0.25, P = P_values$female),
      hsq_F_L20 = mapToLiabilityScale(hsq_female, K = 0.20, P = P_values$female),
      hsq_F_L15 = mapToLiabilityScale(hsq_female, K = 0.15, P = P_values$female),
      hsq_F_L10 = mapToLiabilityScale(hsq_female, K = 0.10, P = P_values$female),
      hsq_F_L5 = mapToLiabilityScale(hsq_female, K = 0.05, P = P_values$female),
      hsq_M_L25 = mapToLiabilityScale(hsq_male, K = 0.25, P = P_values$male),
      hsq_M_L20 = mapToLiabilityScale(hsq_male, K = 0.20, P = P_values$male),
      hsq_M_L15 = mapToLiabilityScale(hsq_male, K = 0.15, P = P_values$male),
      hsq_M_L10 = mapToLiabilityScale(hsq_male, K = 0.10, P = P_values$male),
      hsq_M_L5 = mapToLiabilityScale(hsq_male, K = 0.05, P = P_values$male)
    )
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
AGDS_all_mcmc <- calculate_liability(AGDS_all_mcmc, "AGDS", P_values$AGDS)
AllOfUs_all_mcmc <- calculate_liability(AllOfUs_all_mcmc, "AllOfUs", P_values$AllOfUs)
Bionic_all_mcmc <- calculate_liability(Bionic_all_mcmc, "Bionic", P_values$Bionic)
Blokland21_all_mcmc <- calculate_liability(Blokland21_all_mcmc, "Blokland21", P_values$Blokland21)
GLAD_all_mcmc <- calculate_liability(GLAD_all_mcmc, "GLAD", P_values$GLAD)
UKB_all_mcmc <- calculate_liability(UKB_all_mcmc, "UKB", P_values$UKB)


# Create one df with data from all cohorts
Cohorts_all_mcmc <- bind_rows(
  AGDS_all_mcmc,
  AllOfUs_all_mcmc,
  Bionic_all_mcmc,
  Blokland21_all_mcmc,
  GLAD_all_mcmc,
  UKB_all_mcmc
)




# Create a list of column names, sex labels, and K values for each column
liability_columns <- c("hsq_F_L25", "hsq_F_L20", "hsq_F_L15", "hsq_F_L10", "hsq_F_L5",
                       "hsq_M_L25", "hsq_M_L20", "hsq_M_L15", "hsq_M_L10", "hsq_M_L5")

sex_labels <- c(rep("female", 5), rep("male", 5))
K_values <- c(0.25, 0.20, 0.15, 0.10, 0.05, 0.25, 0.20, 0.15, 0.10, 0.05)




# Initialize an empty list to store results for liability scale
results_list <- list()

# Iterate over each column and calculate posterior mean and 95% HPD interval for liability scale
for (i in seq_along(liability_columns)) {
  column_name <- liability_columns[i]
  sex <- sex_labels[i]
  K <- K_values[i]
  
  # Calculate the posterior mean for liability scale
  # 2,500 data points = chain length of 25,000 with thinning of 10
  # burnin = 5,000 (so 500 data point when thinning of 10)
  posterior_mean_liability <- Cohorts_all_mcmc %>%
    summarize(mean = mean(.data[[column_name]])) %>%
    pull(mean)
  
  # Convert to MCMC object and calculate 95% HPD interval for liability scale
  mcmc_obj_liability <- as.mcmc(Cohorts_all_mcmc %>% pull(.data[[column_name]]))
  hpd_interval_liability <- HPDinterval(mcmc_obj_liability, prob = 0.95)[1, ]
  
  # Store liability scale results in the list
  results_list[[i]] <- data.frame(
    mean_hsq = posterior_mean_liability,
    hpdi_lower = hpd_interval_liability[1],
    hpdi_upper = hpd_interval_liability[2],
    sex = sex,
    K = K
  )
}

# Combine all liability-scale results into a single data frame
h2_liability_varying_K_cohorts <- do.call(rbind, results_list)







### Observed-scale heritability for female and male ###
# Female observed scale
observed_female <- Cohorts_all_mcmc %>%
  summarize(mean = mean(hsq_female)) %>%
  pull(mean)

observed_female_hpd <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(hsq_female)), prob = 0.95)[1, ]

# Male observed scale
observed_male <- Cohorts_all_mcmc %>%
  summarize(mean = mean(hsq_male)) %>%
  pull(mean)

observed_male_hpd <- HPDinterval(as.mcmc(Cohorts_all_mcmc %>% pull(hsq_male)), prob = 0.95)[1, ]

# Create data frame for observed values and bind to liability scale results
observed_h2_cohorts <- data.frame(
  mean_hsq = c(observed_female, observed_male),
  hpdi_lower = c(observed_female_hpd[1], observed_male_hpd[1]),
  hpdi_upper = c(observed_female_hpd[2], observed_male_hpd[2]),
  sex = c("female", "male"),
  K = "observed"
)


# Save the combined data frame
save(h2_liability_varying_K_cohorts, file = paste0(directory, "Cohorts_h2_liability_varying_K.RData"))
save(observed_h2_cohorts, file = paste0(directory, "Cohorts_observed_h2.RData"))

write.table(h2_liability_varying_K_cohorts, file = paste0(directory, "Cohorts_h2_liability_varying_K.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(observed_h2_cohorts, file = paste0(directory, "Cohorts_observed_h2.txt"), col.names = T, row.names = F, sep = "\t", quote = F)



### Forest Plot liability h2 with varying K ###

h2_varying_K_plot_cohorts <- ggplot((h2_liability_varying_K_cohorts %>% filter(K != 0.0563 & K != 0.1018)), aes(x = K, y = mean_hsq*100, colour = sex)) +
  geom_vline(xintercept = 0.10,
             colour = "grey30",
             linetype = 'dashed') +
  geom_vline(xintercept = 0.20,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 1, end = 0) +
  scale_x_continuous("Population Prevalence",
                     breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     labels = c(0.05, 0.10, 0.15, 0.20, 0.25)) +
  scale_y_continuous(expression(italic(h)^2 * " (%)")) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_K_plot_cohorts


ggsave(h2_varying_K_plot_cohorts, file = paste0(directory, "h2_liability_M_F_varying_K_plot_SBayesS_Cohorts.png", sep = ""), width = 10, height = 10, unit = "cm")




### Male vs Female: Using all cohorts ###
# posterior prob that h2 female is > h2 male

# Count the frequency in the MCMC sample that h2_female > h2_male = posterior probability that female h2 is larger than male h2
post_prob_M_F <- Cohorts_all_mcmc %>%
  mutate(hsq_f_minus_m = hsq_female - hsq_male) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

post_prob_M_F_liability <- Cohorts_all_mcmc %>%
  mutate(hsq_f_minus_m = hsq_F_L20 - hsq_M_L10) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_h2_Cohorts.txt"))
cat("Posterior probability that h2_female (observed scale) > h2_male (observed scale):", post_prob_M_F, "\n")
cat("Posterior probability that h2_female (liability scale with K = 0.20) > h2_male (liability scale with K = 0.10):", post_prob_M_F_liability, "\n")
sink()


# Save UK B MCMC samples for h2 so can use in sensitivity analysis of full UKB and downsampled UKB
save(UKB_all_mcmc, file = paste0(directory, "UKB_all_MCMC_female_male.RData"))

write.table(UKB_all_mcmc, file = paste0(directory, "UKB_all_MCMC_female_male.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
