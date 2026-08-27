directory = "/path/10_Publication_figures/03_GxS_Manhattan_plot/"

# GxS = summary statistics file available from download from GWASCatalog
GxS="/path/03_Metal/GxS/fulldc/Metaanalysis_AD_GxS_fulldc_AllCohorts_QCed_rsID.txt"

clumpfile_inter = "/path/05_Clumping/GxS/fulldc/clumped_GxS_fulldc_all_results.txt"


### Packages ###
library(patchwork)
library(showtext)
showtext_auto()
library(tidyverse)

### Load Data ###
GxS_df <- read.table(GxS, header = TRUE, stringsAsFactors = F)

clump_df_inter <- read.table(clumpfile_inter, header = TRUE, stringsAsFactors = F)


### Manhattan Plot ###

# Create cumulative BP position for x axis
cum_BP <- GxS_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

Man_plot_df <- GxS_df %>%
  inner_join(cum_BP, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

# Set significance threshold
sig <- 5e-08
sig2 <- 1e-06


# create df with independent clumps of nominally significant SNPs
clump_sig <- clump_df_inter %>%
  filter(P < 1e-06) %>%
  mutate(SP2 = str_replace_all(SP2, "\\(1\\)", ""),
         SP2 = na_if(SP2, "NONE"))

clump_sig_long <- clump_sig %>%
  separate_longer_delim(SP2, delim = ",")

# create one col of all SNPs
clumps1 <- clump_sig %>%
  select(SNP)

clumps2 <- clump_sig_long %>%
  select(SP2)

clumps_sig_SNPs <- clumps1 %>%
  full_join(clumps2, by = join_by(SNP == SP2)) %>%
  mutate(SNP_in_sig_clump = "TRUE")

# Now in manhattan plot df, create new column with TRUE is SNP is in the list of SNPs part of a significant clump
Man_plot_df <- Man_plot_df %>%
  full_join(clumps_sig_SNPs, by = join_by(MarkerName == SNP)) %>%
  drop_na(-SNP_in_sig_clump)

# Add new column with title of plot so makes title on right like done for Miami plot with facet
Man_plot_df <- Man_plot_df %>%
  mutate(Dataset = "Genotype-by-Sex Interaction")

## Man Plot: P-value GxE ##

# Set max y axis value
ylimP <- Man_plot_df %>%
  filter(P == min(P)) %>%
  slice(1) %>%
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)


# Use subset of SNPs to trial first
Man_plot_df_subset <- Man_plot_df %>%
filter(P < 5e-05)

# Use annotate column to label the lead nominally sig SNPs
clumpinter_sig_list <- clump_sig %>%
  pull(SNP)

Man_plot_df_subset <- Man_plot_df_subset %>%
  mutate(annotate = case_when(
    MarkerName %in% clumpinter_sig_list ~ rsID_build37
  ))


# Plot
Man_plot <- ggplot(Man_plot_df,
                   aes(x = bp_cum, y = -log10(P), color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig),
             color = "grey30",
             linetype = "dashed") +
  geom_hline(yintercept = -log10(sig2),
             color = "grey80",
             linetype = "dashed") +
  geom_point(size = 0.9, alpha = 0.75) +
  geom_point(data = subset(Man_plot_df, SNP_in_sig_clump == TRUE), color = "#481567FF", size = 0.9) +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, ylimP),
                     breaks = seq(0, ylimP, 2),
                     labels = abs(seq(0, ylimP, 2))) +
  scale_color_manual(values = rep(c("#20A387FF", "#95D840FF"),
                                  unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome",
       y = expression(-log[10](p-value))) +
  geom_text(data = subset(Man_plot_df_subset, !is.na(annotate)),
            aes(label = annotate,
                vjust = case_when(
                  annotate == "rs12092435" ~ 0.5,
                  annotate == "rs12312238" ~ -0.5,
                  annotate == "rs6080675" ~ 1,
                  TRUE ~ 0.5
                ),
                hjust = case_when(
                  annotate == "rs12092435" ~ -0.2,
                  annotate == "rs12312238" ~ 1.2,
                  annotate == "rs6080675" ~ -0.2,
                  TRUE ~ 1.2
                )),
            size = 3, family = "Calibri", color = "#481567FF") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(angle = 60, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12, colour = "black"))


# Man_plot



outfile <- paste(directory, "Manhattan_GxS.png", sep="")
ggsave(Man_plot, width = 21, height = 10, unit = "cm", file = outfile)
