args = commandArgs(trailingOnly=TRUE)
infile = args[1]
infile2 = args[2]
facetlabel1 = args[3]
facetlabel2 = args[4]
gwaspw_model1 = args[5]
gwaspw_model2 = args[6]
gwaspw_model3 = args[7]
gwaspw_model4 = args[8]
p_threshold = as.numeric(args[9])
p_threshold2 = as.numeric(args[10])
filename = args[11]

#########################################################################

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(ggh4x)        # To create different scales on ggplot facets
library(tidyverse)

### Load Data ###
GWAS_df <- read.table(infile, header = T, stringsAsFactors = F)

GWAS_df2 <- read.table(infile2, header = T, stringsAsFactors = F)

# Combine the data from GWAS_df and GWAS_df2 into a single dataframe
combined_df <- bind_rows(
  mutate(GWAS_df, Dataset = "GWAS_df"),
  mutate(GWAS_df2, Dataset = "GWAS_df2")
)

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
  mutate(ylim = abs(floor(log10(P))) + 1) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "GWAS_df2") %>%
  filter(P == min(P)) %>%
  mutate(ylim = (abs(floor(log10(P))) + 1)) %>%
  pull(ylim)

# Use the max ylim so the two plots are comparable
chosen_ylim <- max(ylim, ylim2)

# Set significance threshold
sig <- p_threshold
sig2 <- p_threshold2

# gwas-pw results: model 1
model1 <- read.table(gwaspw_model1, header = T, stringsAsFactors = F)

# Create cumulative BP position for x axis
model1_add <- model1 %>%
  mutate(CHR = as.numeric(str_replace(chr, "chr", ""))) %>% 
  left_join(cum_BP)

# Calculate start and stop locations for each segment based on cumulative BP position
# Add dataset label - model1 is evidence for only in 1st sumstats and not 2nd, so add Dataset GWAS_df so appears only in top facet
model1_df <- model1_add %>%
  mutate(start = st + bp_add,
         stop = sp + bp_add,
         Dataset = "GWAS_df")

# gwas-pw results: model 2
model2 <- read.table(gwaspw_model2, header = T, stringsAsFactors = F)

# Create cumulative BP position for x axis
model2_add <- model2 %>%
  mutate(CHR = as.numeric(str_replace(chr, "chr", ""))) %>% 
  left_join(cum_BP)

# Calculate start and stop locations for each segment based on cumulative BP position
# Add dataset label - model2 is evidence for only in 2nd sumstats and not 1st, so add Dataset GWAS_df2 so appears only in bottom facet
model2_df <- model2_add %>%
  mutate(start = st + bp_add,
         stop = sp + bp_add,
         Dataset = "GWAS_df2")


# gwas-pw results: model 3
model3 <- read.table(gwaspw_model3, header = T, stringsAsFactors = F)

# Create cumulative BP position for x axis
model3_add <- model3 %>%
  mutate(CHR = as.numeric(str_replace(chr, "chr", ""))) %>% 
  left_join(cum_BP)

# Calculate start and stop locations for each segment based on cumulative BP position
# model 3 is evidence for both sumstats so don't add dataset and appears in top and bottom facet
model3_df <- model3_add %>%
  mutate(start = st + bp_add,
         stop = sp + bp_add)

# gwas-pw results: model 4
model4 <- read.table(gwaspw_model4, header = T, stringsAsFactors = F)

# Create cumulative BP position for x axis
model4_add <- model4 %>%
  mutate(CHR = as.numeric(str_replace(chr, "chr", ""))) %>% 
  left_join(cum_BP)

# Calculate start and stop locations for each segment based on cumulative BP position
# model 4 is evidence for two distinct associations, so don't add dataset and appears in top and bottom facet (will use different colour in plot)
model4_df <- model4_add %>%
  mutate(start = st + bp_add,
         stop = sp + bp_add)

# Plot (only include gwas-pw results in plot if there is data in each data frame)
facet_labels <- c(
  GWAS_df = paste(facetlabel1),
  GWAS_df2 = paste(facetlabel2)
)

Miami_plot <- ggplot(Man_plot_df,
                      aes(x = bp_cum,
                          y = ifelse(Dataset == "GWAS_df",
                                     -log10(P), log10(P)),
                          color = as_factor(CHR)))

if (nrow(model1_df) > 0) {
  Miami_plot <- Miami_plot +
    geom_rect(data = model1_df,
              aes(xmin = start, xmax = stop),
              ymin = 0, ymax = chosen_ylim,
              fill = "#FDE725FF", colour = "#FDE725FF", alpha = 0.5,
              inherit.aes = FALSE)
}

if (nrow(model2_df) > 0) {
  Miami_plot <- Miami_plot +
    geom_rect(data = model2_df,
              aes(xmin = start, xmax = stop),
              ymin = -chosen_ylim, ymax = 0,
              fill = "#440154FF", colour = "#440154FF", alpha = 0.5,
              inherit.aes = FALSE)
}

if (nrow(model3_df) > 0) {
  Miami_plot <- Miami_plot +
    geom_rect(data = model3_df,
              aes(xmin = start, xmax = stop),
              ymin = -chosen_ylim, ymax = chosen_ylim,
              fill = "grey", colour = "grey", alpha = 0.5,
              inherit.aes = FALSE)
}

if (nrow(model4_df) > 0) {
  Miami_plot <- Miami_plot +
    geom_rect(data = model4_df,
              aes(xmin = start, xmax = stop),
              ymin = -chosen_ylim, ymax = chosen_ylim,
              fill = "red", colour = "red", alpha = 0.5,
              inherit.aes = FALSE)
}

Miami_plot <- Miami_plot +
  geom_point(size = 0.6, alpha = 0.75) +
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "GWAS_df",
                                     -log10(sig), log10(sig))),
             color = "grey30",
             linetype = "dashed") +
  geom_hline(data = Man_plot_df,
             aes(yintercept = ifelse(Dataset == "GWAS_df",
                                     -log10(sig2), log10(sig2))),
             color = "grey80",
             linetype = "dashed") +
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(-chosen_ylim, chosen_ylim, 1),
                     labels = abs(seq(-chosen_ylim, chosen_ylim, 1))) +
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

# To set y axis scales for the two facets with the same chosen_ylim so have same y axis making it easier to compare
position_scales <- list(
  scale_y_continuous(limits = c(0, chosen_ylim),
                     breaks = seq(0, chosen_ylim, 1),
                     labels = abs(seq(0, chosen_ylim, 1))),
  scale_y_continuous(limits = c(-chosen_ylim, 0),
                     breaks = seq(-chosen_ylim, 0, 1),
                     labels = abs(seq(-chosen_ylim, 0, 1))))

Miami_plot <- Miami_plot + facetted_pos_scales(y = position_scales)

outfile <- paste(filename, ".Miami.gwaspw.png", sep="")
ggsave(Miami_plot, width = 40, height = 20, unit = "cm", file = outfile)

