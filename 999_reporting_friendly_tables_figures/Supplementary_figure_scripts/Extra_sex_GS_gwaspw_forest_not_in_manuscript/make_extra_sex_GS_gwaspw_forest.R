# Create supplementary figure - forest plot of effect sizes from female and male GWAS for SNPs identified as genome-wide significant in the female GWAS, male GWAS, nominally significant in the GxS analysis and possible causal variants from gwas-pw

directory = "/path/10_Publication_figures/Supplementary_figures/Forest_plots_female_male_GxS_gwaspw/"

GxS="/path/03_Metal/GxS/fulldc/Metaanalysis_AD_GxS_fulldc_AllCohorts_QCed_rsID.txt"

clumpfile_inter = "/path/05_Clumping/GxS/fulldc/clumped_GxS_fulldc_all_results.txt"

femaleGWAS = "/path/03_Metal/Females/Metaanalysis_AD_female_AllCohorts_QCed_rsID.txt"

maleGWAS = "/path/03_Metal/Males/Metaanalysis_AD_male_AllCohorts_QCed_rsID.txt"

clumpfile_female = "/path/05_Clumping/Females/clumped_female_all_results.txt"

clumpfile_male = "/path/05_Clumping/Males/clumped_male_all_results.txt"


gwaspw_causal_shared = "/path/06_GWAS_pw/femaleAD_maleAD/SNP_with_PPA_above_05_per_shared_region.txt"

### Packages ###
library(patchwork)
library(tidyverse)

### Load Data ###
GxS_df <- read.table(GxS, header = TRUE, stringsAsFactors = F)

clump_df_inter <- read.table(clumpfile_inter, header = TRUE, stringsAsFactors = F)
                             
femaleGWAS_df <- read.table(femaleGWAS, header = T, stringsAsFactors = F)
                             
maleGWAS_df <- read.table(maleGWAS, header = T, stringsAsFactors = F)
                             
clump_df_female <- read.table(clumpfile_female, header = TRUE, stringsAsFactors = F)
                             
clump_df_male <- read.table(clumpfile_male, header = TRUE, stringsAsFactors = F)

gwaspw_causal_shared_df <- read.table(gwaspw_causal_shared, header = T, stringsAsFactors = F)
                             
                             
### FOREST PLOTS ###
                            
# FOREST PLOT: FEMALE SNPs
# Independent lead genome-wide sig SNPs from female GWAS
clump_female_sig <- clump_df_female %>% 
  filter(P < 5e-08) %>%
  select(SNP) %>%
  rename(MarkerName = SNP)

# Extract these SNPs for both female and male sumstas
femaleGWAS_female_sig <- femaleGWAS_df %>%
  inner_join(clump_female_sig) %>%
  mutate(Dataset = "female")

maleGWAS_female_sig <- maleGWAS_df %>%
  inner_join(clump_female_sig) %>%
  mutate(Dataset = "male")

# Create one df with all data to make plot
df_plot_female <- femaleGWAS_female_sig %>%
  full_join(maleGWAS_female_sig)

# add CI                        
df_plot_female <- df_plot_female %>%
  mutate(MarkerName = factor(MarkerName, levels = unique(MarkerName)),
         CI = 1.96 * StdErr)

# Order SNPs by position
df_plot_female <- df_plot_female %>% 
  mutate(rsID_build37 = factor(rsID_build37,
                               levels = df_plot_female %>%
                                 distinct(rsID_build37, CHR, BP) %>%
                                 arrange(CHR, BP) %>%
                                 pull(rsID_build37) %>% 
                                 rev()))
# Plot                           
forest_plot_female <- ggplot(df_plot_female, aes(x = Effect, y = rsID_build37, colour = Dataset)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dotted') +
  geom_linerange(aes(xmin = Effect - CI, xmax = Effect + CI),
                 position = position_dodge(width = -0.8)) +
  geom_point(position = position_dodge(width = -0.8)) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Females", "Males")) +
  scale_x_continuous("Effect Size",
                     lim = c(-0.10, 0.15),
                     breaks = seq(-0.1, 0.15, 0.05),
                     labels = seq(-0.1, 0.15, 0.05)) +
  scale_y_discrete("") +
  labs(title = "Female SNPs") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(0,10,0,0), "pt"),
        plot.title = element_text(size = 12))

forest_plot_female


# FOREST PLOT: MALE SNPs
# Independent lead genome-wide sig SNPs from male GWAS
clump_male_sig <- clump_df_male %>% 
  filter(P < 5e-08) %>%
  select(SNP) %>%
  rename(MarkerName = SNP)

# Extract these SNPs for both female and male sumstas
femaleGWAS_male_sig <- femaleGWAS_df %>%
  inner_join(clump_male_sig) %>%
  mutate(Dataset = "female")

maleGWAS_male_sig <- maleGWAS_df %>%
  inner_join(clump_male_sig) %>%
  mutate(Dataset = "male")

# Create one df with all data to make plot
df_plot_male <- femaleGWAS_male_sig %>%
  full_join(maleGWAS_male_sig)

# add CI                       
df_plot_male <- df_plot_male %>%
  mutate(MarkerName = factor(MarkerName, levels = unique(MarkerName)),
         CI = 1.96 * StdErr)

# Order SNPs by position
df_plot_male <- df_plot_male %>% 
  mutate(rsID_build37 = factor(rsID_build37,
                               levels = df_plot_male %>%
                                 distinct(rsID_build37, CHR, BP) %>%
                                 arrange(CHR, BP) %>%
                                 pull(rsID_build37) %>% 
                                 rev()))
# Plot                           
forest_plot_male <- ggplot(df_plot_male, aes(x = Effect, y = rsID_build37, colour = Dataset)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dotted') +
  geom_linerange(aes(xmin = Effect - CI, xmax = Effect + CI),
                 position = position_dodge(width = -0.8)) +
  geom_point(position = position_dodge(width = -0.8)) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Females", "Males")) +
  scale_x_continuous("Effect Size",
                     lim = c(-0.10, 0.15),
                     breaks = seq(-0.1, 0.15, 0.05),
                     labels = seq(-0.1, 0.15, 0.05)) +
  scale_y_discrete("") +
  labs(title = "Male SNPs") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(0,10,0,0), "pt"),
        plot.title = element_text(size = 12))

forest_plot_male



# FOREST PLOT: GxS SNPs 
# Independent lead nominally sig SNPs from GxS analysis
clumpinter_sig <- clump_df_inter %>%
  filter(P < 1e-06) %>%
  select(SNP) %>%
  rename(MarkerName = SNP)

# Extract these SNPs for both female and male sumstas
femaleGWAS_clumpinter_sig <- femaleGWAS_df %>%
  inner_join(clumpinter_sig) %>%
  mutate(Dataset = "female")

maleGWAS_df_clumpinter_sig <- maleGWAS_df %>%
  inner_join(clumpinter_sig) %>%
  mutate(Dataset = "male")

# Create one df with all data to make plot
df_plot_GxS <- femaleGWAS_clumpinter_sig %>%
  full_join(maleGWAS_df_clumpinter_sig)

# add CI                      
df_plot_GxS <- df_plot_GxS %>%
  mutate(MarkerName = factor(MarkerName, levels = unique(MarkerName)),
         CI = 1.96 * StdErr)

# Order SNPs by position
df_plot_GxS <- df_plot_GxS %>% 
  mutate(rsID_build37 = factor(rsID_build37,
                               levels = df_plot_GxS %>%
                                 distinct(rsID_build37, CHR, BP) %>%
                                 arrange(CHR, BP) %>%
                                 pull(rsID_build37) %>% 
                                 rev()))
                             
# Plot                           
forest_plot_GxS <- ggplot(df_plot_GxS, aes(x = Effect, y = rsID_build37, colour = Dataset)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dotted') +
  geom_linerange(aes(xmin = Effect - CI, xmax = Effect + CI),
                 position = position_dodge(width = -0.8)) +
  geom_point(position = position_dodge(width = -0.8)) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Females", "Males")) +
  scale_x_continuous("Effect Size",
                     lim = c(-0.10, 0.1),
                     breaks = seq(-0.1, 0.1, 0.05),
                     labels = seq(-0.1, 0.1, 0.05)) +
  scale_y_discrete("") +
  labs(title = "GxS SNPs") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(0,10,0,0), "pt"),
        plot.title = element_text(size = 12))
                             
forest_plot_GxS



# FOREST PLOT: GWASPW POSSIBLE CAUSAL VARIANTS IN BOTH SEXES
# SNPs identified in gwas-pw as possible causal varints for AD in both sexes
gwaspw_causal_shared_df <- gwaspw_causal_shared_df %>% 
  rename(MarkerName = id)

# Extract these SNPs for both female and male sumstas
femaleGWAS_gwaspw <- femaleGWAS_df %>%
  inner_join(gwaspw_causal_shared_df) %>%
  mutate(Dataset = "female")

maleGWAS_gwaspw <- maleGWAS_df %>%
  inner_join(gwaspw_causal_shared_df) %>%
  mutate(Dataset = "male")

# Create one df with all data to make plot
df_plot_gwaspw <- femaleGWAS_gwaspw %>%
  full_join(maleGWAS_gwaspw)

# add CI                      
df_plot_gwaspw <- df_plot_gwaspw %>%
  mutate(MarkerName = factor(MarkerName, levels = unique(MarkerName)),
         CI = 1.96 * StdErr)

# Order SNPs by position
df_plot_gwaspw <- df_plot_gwaspw %>% 
  mutate(rsID_build37 = factor(rsID_build37,
                               levels = df_plot_gwaspw %>%
                                 distinct(rsID_build37, CHR, BP) %>%
                                 arrange(CHR, BP) %>%
                                 pull(rsID_build37) %>% 
                                 rev()))
# Plot                           
forest_plot_gwaspw <- ggplot(df_plot_gwaspw, aes(x = Effect, y = rsID_build37, colour = Dataset)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dotted') +
  geom_linerange(aes(xmin = Effect - CI, xmax = Effect + CI),
                 position = position_dodge(width = -0.8)) +
  geom_point(position = position_dodge(width = -0.8)) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Females", "Males")) +
  scale_x_continuous("Effect Size",
                     lim = c(-0.10, 0.15),
                     breaks = seq(-0.1, 0.15, 0.05),
                     labels = seq(-0.1, 0.15, 0.05)) +
  scale_y_discrete("") +
  labs(title = "Possible causal variants in both sexes (gwas-pw)") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(0,10,0,0), "pt"),
        plot.title = element_text(size = 12))

forest_plot_gwaspw


# COMBINED FIGURE
figure <- (forest_plot_female | forest_plot_male) /
  (forest_plot_GxS | forest_plot_gwaspw) +
  plot_layout(heights = c(0.7, 0.3)) +
  plot_annotation(tag_levels = 'A')

figure
                           
outfile <- paste(directory, "Forest_plots.png", sep="")
ggsave(figure, width = 22, height = 22, unit = "cm", file = outfile)
                             
