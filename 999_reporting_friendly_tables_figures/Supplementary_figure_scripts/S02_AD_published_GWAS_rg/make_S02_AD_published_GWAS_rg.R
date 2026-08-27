directory = "/path/10_Publication_figures/Supplementary_figures/Rg_females_males_AD/"

rg_female="/path/04_LDSC/SNPrg/Females/rg_femaleAD_correlation_table.txt"

rg_male="/path/04_LDSC/SNPrg/Males/rg_maleAD_correlation_table.txt"

### LOAD PACKAGES ###

library(patchwork)
library(tidyverse)



### LOAD DATA ###

rg_female_df <- read.table(rg_female, header = T, stringsAsFactors = F)

rg_male_df <- read.table(rg_male, header = T, stringsAsFactors = F)

rg_df <- bind_rows(rg_female_df, rg_male_df) %>%
  mutate(CI = 1.96 * se,
         p1 = gsub("_sumstats", "", p1),
         p2 = gsub("_sumstats", "", p2))


### RG WITH PREVIOUS AD ###

order <- c('AdamsAD',
           'BloklandADfemale',
           'BloklandADmale',
           'SilveiraADfemale',
           'SilveiraADmale')
order <- rev(order)

plot_df <- rg_df %>%
  filter(p2 == 'AdamsAD' |
           p2 == 'BloklandADfemale' |
           p2 == 'BloklandADmale' |
           p2 == 'SilveiraADfemale' |
           p2 == 'SilveiraADmale') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_AD", "Allcohorts_male_AD")))



rg_plot_AD <- ggplot(plot_df, aes(x = rg, y = p2, colour = p1)) +
  geom_pointrange(aes(x = rg, xmin = rg + CI, xmax = rg - CI),
                  position = position_dodge(width = -0.8)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete(breaks = c('AdamsAD',
                              'BloklandADfemale',
                              'BloklandADmale',
                              'SilveiraADfemale',
                              'SilveiraADmale'),
                   labels = c('Adams et al., 2025',
                              'Blokland et al., 2022: Females',
                              'Blokland et al., 2022: Males',
                              'Silveira et al., 2023: Females',
                              'Silveira et al., 2023: Males')) +
  labs(x = "Genetic Correlation",
       y = "") +
  annotate("text", x = 1.1, y = "AdamsAD", label = "*", size = 5) +
  annotate("text", x = 1.3, y = "BloklandADfemale", label = "*", size = 5) +
  annotate("text", x = 1.7, y = "BloklandADmale", label = "*", size = 5) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        legend.position = "top",
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))


rg_plot_AD

outfile <- paste(directory, "Rg_female_vs_male_AD.png", sep="")
ggsave(rg_plot_AD, width = 20, height = 10, unit = "cm", file = outfile)
