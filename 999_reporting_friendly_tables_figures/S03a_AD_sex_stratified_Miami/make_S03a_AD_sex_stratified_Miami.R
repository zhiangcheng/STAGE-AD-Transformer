directory = "/path/10_Publication_figures/01_Sex_stratified_Miami_plot/"

# femaleGWAS and maleGWAS are the sex-stratified summary statistics files available to download from GWASCatalog
femaleGWAS = "/path/03_Metal/Females/Metaanalysis_AD_female_AllCohorts_QCed_rsID.txt"

maleGWAS = "/path/03_Metal/Males/Metaanalysis_AD_male_AllCohorts_QCed_rsID.txt"

clumpfile_female = "/path/05_Clumping/Females/clumped_female_all_results.txt"

clumpfile_male = "/path/05_Clumping/Males/clumped_male_all_results.txt"


#########################################################################

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(ggh4x)        # To create different scales on ggplot facets
library(showtext)
showtext_auto()
library(tidyverse)

### Load Data ###
femaleGWAS_df <- read.table(femaleGWAS, header = T, stringsAsFactors = F)

maleGWAS_df <- read.table(maleGWAS, header = T, stringsAsFactors = F)

clump_df_female <- read.table(clumpfile_female, header = TRUE, stringsAsFactors = F)

clump_df_male <- read.table(clumpfile_male, header = TRUE, stringsAsFactors = F)

# For females create df with independent clumps of genome-wide significant SNPs
clump_sig_female <- clump_df_female %>% 
  filter(P < 5e-08) %>% 
  mutate(SP2 = str_replace_all(SP2, "\\(1\\)", ""),
         SP2 = na_if(SP2, "NONE"))

clump_sig_female_long <- clump_sig_female %>% 
  separate_longer_delim(SP2, delim = ",")

# create one col of all SNPs
clumps1_female <- clump_sig_female %>% 
  select(SNP)

clumps2_female <- clump_sig_female_long %>% 
  select(SP2)

clumps_sig_SNPs_female <- clumps1_female %>% 
  full_join(clumps2_female, by = join_by(SNP == SP2)) %>% 
  mutate(SNP_in_sig_clump = "TRUE")

# Now in GWAS df, create new column with TRUE is SNP is in the list of SNPs part of a significant clump
femaleGWAS_df_clumps <- femaleGWAS_df %>% 
  full_join(clumps_sig_SNPs_female, by = join_by(MarkerName == SNP)) %>% 
  drop_na(-SNP_in_sig_clump) 

# For males create df with independent clumps of genome-wide significant SNPs
clump_sig_male <- clump_df_male %>% 
  filter(P < 5e-08) %>% 
  mutate(SP2 = str_replace_all(SP2, "\\(1\\)", ""),
         SP2 = na_if(SP2, "NONE"))

clump_sig_male_long <- clump_sig_male %>% 
  separate_longer_delim(SP2, delim = ",")

# create one col of all SNPs
clumps1_male <- clump_sig_male %>% 
  select(SNP)

clumps2_male <- clump_sig_male_long %>% 
  select(SP2)

clumps_sig_SNPs_male <- clumps1_male %>% 
  full_join(clumps2_male, by = join_by(SNP == SP2)) %>% 
  mutate(SNP_in_sig_clump = "TRUE")

# Now in GWAS df, create new column with TRUE is SNP is in the list of SNPs part of a significant clump
maleGWAS_df_clumps <- maleGWAS_df %>% 
  full_join(clumps_sig_SNPs_male, by = join_by(MarkerName == SNP)) %>% 
  drop_na(-SNP_in_sig_clump) 



# Combine the data from females and males into a single dataframe
combined_df <- bind_rows(
  mutate(femaleGWAS_df_clumps, Dataset = "female"),
  mutate(maleGWAS_df_clumps, Dataset = "male")
)

# Filter data below a certain p-value so runs fast as a quick check
# combined_df <- combined_df %>%
#   filter(P < 5e-05)

# Create cumulative BP position for x axis
cum_BP <- combined_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

Man_plot_df <- combined_df %>%
  inner_join(cum_BP, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

# Set max y axis value
ylim <- Man_plot_df %>%
  filter(Dataset == "female") %>%
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "male") %>%
  filter(P == min(P)) %>%
  mutate(ylim = (abs(floor(log10(P))) + 1)) %>%
  pull(ylim)

# Use the max ylim so the two plots are comparable
chosen_ylim <- max(ylim, ylim2)

# Set significance threshold
sig <- 5e-08
sig2 <- 1e-06

# Plot
facet_labels <- c(
  female = paste("Females"),
  male = paste("Males")
)


Miami_plot <- ggplot(Man_plot_df,
                     aes(x = bp_cum,
                         y = ifelse(Dataset == "female",
                                    -log10(P), log10(P)),
                         color = as_factor(CHR)))+
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "female",
                                     -log10(sig), log10(sig))),
             color = "grey30",
             linetype = "dashed") +
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "female",
                                     -log10(sig2), log10(sig2))),
             color = "grey80",
             linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75) +
  geom_point(data = subset(Man_plot_df, SNP_in_sig_clump == TRUE), color = "#481567FF", size = 0.9, alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(-chosen_ylim, chosen_ylim, 1),
                     labels = abs(seq(-chosen_ylim, chosen_ylim, 1))) +
  scale_color_manual(values = rep(c("#20A387FF", "#95D840FF"),
                                  unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome",
       y = expression(-log[10](paste("p-value"))),
       color = "Dataset") +
  facet_wrap(~ Dataset, scale = "free_y", ncol = 1, strip.position = "right",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, colour = "#191d1f", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 7, colour = "#191d1f", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 7, colour = "#191d1f"),
        axis.line = element_line(colour = "#191d1f"),
        axis.ticks = element_line(colour = "#191d1f"))

# To set y axis scales for the two facets with the same chosen_ylim so have same y axis making it easier to compare
position_scales <- list(
  scale_y_continuous(limits = c(0, chosen_ylim),
                     breaks = seq(0, chosen_ylim, 2),
                     labels = abs(seq(0, chosen_ylim, 2))),
  scale_y_continuous(limits = c(-chosen_ylim, 0),
                     breaks = seq(-chosen_ylim, 0, 2),
                     labels = abs(seq(-chosen_ylim, 0, 2))))

Miami_plot <- Miami_plot + facetted_pos_scales(y = position_scales)

# Miami_plot

# outfile <- paste(directory, "Miami_plot_sex_stratified.png", sep="")
# ggsave(Miami_plot, width = 18, height = 9, unit = "cm", file = outfile)

outfile <- paste(directory, "Miami_plot_sex_stratified.eps", sep="")
ggsave(Miami_plot, width = 18, height = 9, units = "cm", file = outfile, device = cairo_ps, dpi = 300)

