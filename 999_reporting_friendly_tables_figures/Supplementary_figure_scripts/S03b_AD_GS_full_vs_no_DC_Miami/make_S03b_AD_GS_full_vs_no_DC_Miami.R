# Create Miami plot of GxS usung full dc and GxS using no dc

directory = "/path/10_Publication_figures/Supplementary_figures/Miami_GxS_fulldc_nodc/"

fulldcGWAS = "/path/Meta-analysis/Freeze2/03_Metal/GxS/fulldc/Metaanalysis_AD_GxS_fulldc_AllCohorts_QCed_rsID.txt"

nodcGWAS = "/path/Meta-analysis/Freeze2/03_Metal/GxS/nodc/Metaanalysis_AD_GxS_nodc_AllCohorts_QCed_rsID.txt"

clumpfile_fulldc = "/path/Meta-analysis/Freeze2/05_Clumping/GxS/fulldc/clumped_GxS_fulldc_all_results.txt"

clumpfile_nodc = "/path/Meta-analysis/Freeze2/05_Clumping/GxS/nodc/clumped_GxS_nodc_all_results.txt"


#########################################################################

### Packages ###
library(ggh4x)        # To create different scales on ggplot facets
library(tidyverse)

### Load Data ###
fulldcGWAS_df <- read.table(fulldcGWAS, header = T, stringsAsFactors = F)

nodcGWAS_df <- read.table(nodcGWAS, header = T, stringsAsFactors = F)

clump_df_fulldc <- read.table(clumpfile_fulldc, header = TRUE, stringsAsFactors = F)

clump_df_nodc <- read.table(clumpfile_nodc, header = TRUE, stringsAsFactors = F)

# For fulldcs create df with independent clumps of genome-wide significant SNPs
clump_sig_fulldc <- clump_df_fulldc %>% 
  filter(P < 1e-06) %>% 
  mutate(SP2 = str_replace_all(SP2, "\\(1\\)", ""),
         SP2 = na_if(SP2, "NONE"))

clump_sig_fulldc_long <- clump_sig_fulldc %>% 
  separate_longer_delim(SP2, delim = ",")

# create one col of all SNPs
clumps1_fulldc <- clump_sig_fulldc %>% 
  select(SNP)

clumps2_fulldc <- clump_sig_fulldc_long %>% 
  select(SP2)

clumps_sig_SNPs_fulldc <- clumps1_fulldc %>% 
  full_join(clumps2_fulldc, by = join_by(SNP == SP2)) %>% 
  mutate(SNP_in_sig_clump = "TRUE")

# Now in GWAS df, create new column with TRUE is SNP is in the list of SNPs part of a significant clump
fulldcGWAS_df_clumps <- fulldcGWAS_df %>% 
  full_join(clumps_sig_SNPs_fulldc, by = join_by(MarkerName == SNP)) %>% 
  drop_na(-SNP_in_sig_clump) 

# For nodcs create df with independent clumps of genome-wide significant SNPs
clump_sig_nodc <- clump_df_nodc %>% 
  filter(P < 1e-06) %>% 
  mutate(SP2 = str_replace_all(SP2, "\\(1\\)", ""),
         SP2 = na_if(SP2, "NONE"))

clump_sig_nodc_long <- clump_sig_nodc %>% 
  separate_longer_delim(SP2, delim = ",")

# create one col of all SNPs
clumps1_nodc <- clump_sig_nodc %>% 
  select(SNP)

clumps2_nodc <- clump_sig_nodc_long %>% 
  select(SP2)

clumps_sig_SNPs_nodc <- clumps1_nodc %>% 
  full_join(clumps2_nodc, by = join_by(SNP == SP2)) %>% 
  mutate(SNP_in_sig_clump = "TRUE")

# Now in GWAS df, create new column with TRUE is SNP is in the list of SNPs part of a significant clump
nodcGWAS_df_clumps <- nodcGWAS_df %>% 
  full_join(clumps_sig_SNPs_nodc, by = join_by(MarkerName == SNP)) %>% 
  drop_na(-SNP_in_sig_clump) 



# Combine the data from fulldcs and nodcs into a single dataframe
combined_df <- bind_rows(
  mutate(fulldcGWAS_df_clumps, Dataset = "fulldc"),
  mutate(nodcGWAS_df_clumps, Dataset = "nodc")
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
  filter(Dataset == "fulldc") %>%
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "nodc") %>%
  filter(P == min(P)) %>%
  mutate(ylim = (abs(floor(log10(P))) + 1)) %>%
  pull(ylim)

# Use the max ylim so the two plots are comparable
# chosen_ylim <- max(ylim, ylim2)
chosen_ylim <- 8

# Set significance threshold
sig <- 5e-08
sig2 <- 1e-06

# Plot
facet_labels <- c(
  fulldc = paste("Full Dosage Compensation"),
  nodc = paste("No Dosage Compensation")
)


Miami_plot <- ggplot(Man_plot_df,
                     aes(x = bp_cum,
                         y = ifelse(Dataset == "fulldc",
                                    -log10(P), log10(P)),
                         color = as_factor(CHR)))+
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "fulldc",
                                     -log10(sig), log10(sig))),
             color = "grey30",
             linetype = "dashed") +
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "fulldc",
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
       y = expression(-log[10](p-value)),
       color = "Dataset") +
  facet_wrap(~ Dataset, scale = "free_y", ncol = 1, strip.position = "right",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "#191d1f", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "#191d1f", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 10, colour = "#191d1f"),
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

outfile <- paste(directory, "Miami_plot_GxS_fulldc_nodc.png", sep="")
ggsave(Miami_plot, width = 19, height = 10, unit = "cm", file = outfile)


