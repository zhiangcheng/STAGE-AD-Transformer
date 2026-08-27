rg_female="/path/04_LDSC/SNPrg/Females/rg_femaleMDD_correlation_table.txt"

rg_male="/path/04_LDSC/SNPrg/Males/rg_maleMDD_correlation_table.txt"



### LOAD PACKAGES ###

library(tidyverse)



### READ IN DATA ###

rg_female_df <- read.table(rg_female, header = T, stringsAsFactors = F)

rg_male_df <- read.table(rg_male, header = T, stringsAsFactors = F)


### DATAFRAME

rg_df <- bind_rows(rg_female_df, rg_male_df) %>%
  mutate(CI = 1.96 * se,
         p1 = gsub("_sumstats", "", p1),
         p2 = gsub("_sumstats", "", p2))



### RG WITH OTHER TRAITS ###

order <- c('Schizophrenia',
           'PTSD',
           'Bipolar',
           'Anxiety',
           'ADHD',
           'EA',
           'BMI',
           'Waist_to_hip_ratio',
           'Metabolic_syndrome',
           'Drinks_per_week',
           'Ever_smoked_regularly')
order <- rev(order)

plot_df <- rg_df %>%
  filter(p2 == 'Schizophrenia' |
           p2 == 'PTSD' |
           p2 == 'Bipolar' |
           p2 == 'Anxiety' |
           p2 == 'ADHD' |
           p2 == 'EA' |
           p2 == 'BMI' |
           p2 == 'Waist_to_hip_ratio' |
           p2 == 'Metabolic_syndrome' |
           p2 == 'Drinks_per_week' |
           p2 == 'Ever_smoked_regularly') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_MDD", "Allcohorts_male_MDD")))



rg_plot_traits <- ggplot(plot_df, aes(x = rg, y = p2, colour = p1)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_pointrange(aes(x = rg, xmin = rg + CI, xmax = rg - CI),
                  position = position_dodge(width = -0.8)) +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete(breaks = c('Schizophrenia',
                              'PTSD',
                              'Bipolar',
                              'Anxiety',
                              'ADHD',
                              'EA',
                              'BMI',
                              'Waist_to_hip_ratio',
                              'Metabolic_syndrome',
                              'Drinks_per_week',
                              'Ever_smoked_regularly'),
                   labels = c('Schizophrenia',
                              'PTSD',
                              'Bipolar',
                              'Anxiety',
                              'ADHD',
                              'Educational Attainment',
                              'BMI',
                              'Waist/hip ratio',
                              'Metabolic syndrome',
                              'Drinks/week',
                              'Regular smoker')) +
  labs(x = "Genetic Correlation",
       y = "") +
  annotate("text", x = 0.7, y = "ADHD", label = "*", size = 5) +
  annotate("text", x = 0.3, y = "BMI", label = "*", size = 5) +
  annotate("text", x = 0.35, y = "Metabolic_syndrome", label = "*", size = 5) +
  annotate("text", x = 0.45, y = "Ever_smoked_regularly", label = "*", size = 5) +
  theme_classic()

ggsave(rg_plot_traits, width = 20, height = 10, unit = "cm", file = "/path/04_LDSC/SNPrg/forest_plot_rg_traits.png")


### RG WITH SEX-SPECIFIC BMI ###

order <- c('BMI_female',
           'BMI_male')
order <- rev(order)

plot_df <- rg_df %>%
  filter(p2 == 'BMI_female' |
           p2 == 'BMI_male') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_MDD", "Allcohorts_male_MDD")))



rg_plot_BMI <- ggplot(plot_df, aes(x = rg, y = p2, colour = p1)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_pointrange(aes(x = rg, xmin = rg + CI, xmax = rg - CI),
                  position = position_dodge(width = -0.8)) +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete(breaks = c('BMI_female',
                              'BMI_male'),
                   labels = c('BMI: Females',
                              'BMI: Males')) +
  labs(x = "Genetic Correlation",
       y = "") +
  annotate("text", x = 0.4, y = "BMI_female", label = "*", size = 5) +
  annotate("text", x = 0.4, y = "BMI_male", label = "*", size = 5) +
  theme_classic()

ggsave(rg_plot_BMI, width = 20, height = 10, unit = "cm", file = "/path/04_LDSC/SNPrg/forest_plot_rg_BMI_sex_specific.png")


### RG WITH PREVIOUS MDD ###

order <- c('AdamsMDD',
           'BloklandMDDfemale',
           'BloklandMDDmale',
           'SilveiraMDDfemale',
           'SilveiraMDDmale')
order <- rev(order)

plot_df <- rg_df %>%
  filter(p2 == 'AdamsMDD' |
           p2 == 'BloklandMDDfemale' |
           p2 == 'BloklandMDDmale' |
           p2 == 'SilveiraMDDfemale' |
           p2 == 'SilveiraMDDmale') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_MDD", "Allcohorts_male_MDD")))



rg_plot_MDD <- ggplot(plot_df, aes(x = rg, y = p2, colour = p1)) +
  geom_pointrange(aes(x = rg, xmin = rg + CI, xmax = rg - CI),
                  position = position_dodge(width = -0.8)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete(breaks = c('AdamsMDD',
                              'BloklandMDDfemale',
                              'BloklandMDDmale',
                              'SilveiraMDDfemale',
                              'SilveiraMDDmale'),
                   labels = c('Adams et al., 2024',
                              'Blokland et al., 2022: Females',
                              'Blokland et al., 2022: Males',
                              'Silveira et al., 2023: Females',
                              'Silveira et al., 2023: Males')) +
  labs(x = "Genetic Correlation",
       y = "") +
  annotate("text", x = 1.1, y = "AdamsMDD", label = "*", size = 5) +
  annotate("text", x = 1.3, y = "BloklandMDDfemale", label = "*", size = 5) +
  annotate("text", x = 1.7, y = "BloklandMDDmale", label = "*", size = 5) +
  theme_classic()



ggsave(rg_plot_MDD, width = 20, height = 10, unit = "cm", file = "/path/04_LDSC/SNPrg/forest_plot_rg_MDD.png")
