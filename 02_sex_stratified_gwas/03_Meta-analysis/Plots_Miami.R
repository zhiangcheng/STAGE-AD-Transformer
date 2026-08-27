args = commandArgs(trailingOnly=TRUE)
infile = args[1]
infile2 = args[2]
facetlabel1 = args[3]
facetlabel2 = args[4]
filename = args[5]
p_threshold = as.numeric(args[6])

#########################################################################

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(ggh4x)        # To create different scales on ggplot facets
library(tidyverse)

### Load Data ###
GWAS_df <- read.table(infile, header = TRUE, stringsAsFactors = F)

GWAS_df2 <- read.table(infile2, header = TRUE, stringsAsFactors = F)

# Combine the data from GWAS_df and GWAS_df2 into a single dataframe
combined_df <- bind_rows(
  mutate(GWAS_df, Dataset = "GWAS_df"),
  mutate(GWAS_df2, Dataset = "GWAS_df2")
)

### Manhattan Plot ###
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
  filter(Dataset == "GWAS_df") %>% 
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 2) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "GWAS_df2") %>% 
  filter(P == min(P)) %>%
  mutate(ylim = (abs(floor(log10(P))) + 2)) %>%
  pull(ylim)

# Use the max ylim so the two plots are comparable
# Currently the chosen_ylim is not working - plot just uses auto y axis limits. 
# However code in plot scale_y_continuous does make sure y axis is positive values above and below x axis
chosen_ylim <- max(ylim, ylim2)

# Set significance threshold
sig <- p_threshold


# Plot
facet_labels <- c(
  GWAS_df = paste(facetlabel1),
  GWAS_df2 = paste(facetlabel2)
)

(Miami_plot <- ggplot(Man_plot_df,
                      aes(x = bp_cum, 
                          y = ifelse(Dataset == "GWAS_df",
                                     -log10(P), log10(P)), 
                          color = as_factor(CHR))) +
    geom_point(size = 0.6, alpha = 0.75) +
    geom_hline(data = Man_plot_df,
               aes(yintercept = ifelse(Dataset == "GWAS_df",
                                       -log10(sig), log10(sig))),
               color = "grey40",
               linetype = "dashed") +
    scale_x_continuous(label = axis_set$CHR,
                       breaks = axis_set$center) +
    scale_y_continuous(breaks = seq(-chosen_ylim,chosen_ylim,1), 
                       labels = abs(seq(-chosen_ylim,chosen_ylim,1))) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"),
                                    unique(length(axis_set$CHR)))) +
    labs(x = "Chromosome",
         y = expression(log[10](p-value)),
         color = "Dataset") +
    facet_wrap(~ Dataset, scale = "free_y", ncol = 1, strip.position = "right",
               labeller = as_labeller(facet_labels)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
)

# To set y axis scales for the two facets with the same chosen_ylim so have same y axis making it easier to compare
position_scales <- list(
  scale_y_continuous(limits = c(0,chosen_ylim),
                     breaks = seq(-chosen_ylim,chosen_ylim,1),
                     labels = abs(seq(-chosen_ylim,chosen_ylim,1))),
  scale_y_continuous(limits = c(-chosen_ylim,0),
                     breaks = seq(-chosen_ylim,chosen_ylim,1),
                     labels = abs(seq(-chosen_ylim,chosen_ylim,1)))
)

Miami_plot <- Miami_plot + facetted_pos_scales(y = position_scales)

outfile <- paste(filename, ".Miami.png", sep="")
ggsave(Miami_plot, width = 40, height = 20, unit = "cm", file = outfile)


