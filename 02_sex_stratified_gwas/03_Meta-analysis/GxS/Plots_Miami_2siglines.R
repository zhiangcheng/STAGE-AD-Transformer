Sys.setLanguage("en")

args = commandArgs(trailingOnly = TRUE)

infile = "/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/020_Meta_analysis/GxS/fulldc/Metaanalysis_ad_GxS_fulldc_AllCohorts_QCed_rsID.txt"
infile2 = "/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/020_Meta_analysis/GxS/nodc/Metaanalysis_ad_GxS_nodc_AllCohorts_QCed_rsID.txt"

facetlabel1 = "Full_DC"
facetlabel2 = "No_DC"
filename = "fulldc_nodc_QCed_rsID"

p_threshold = 5e-08
p_threshold2 = 1e-06

#########################################################################

### Packages ###
library(data.table)
library(qqman)
library(ggh4x)
library(tidyverse)
library(Rmpfr)

#########################################################################

### Find column name automatically ###

find_col <- function(df, candidates, required_name) {
  
  nms <- colnames(df)
  nms_lower <- tolower(nms)
  candidates_lower <- tolower(candidates)
  
  hit <- which(nms_lower %in% candidates_lower)
  
  if (length(hit) == 0) {
    stop(
      paste0(
        "Cannot find column for ", required_name, ".\n",
        "Tried: ", paste(candidates, collapse = ", "), "\n",
        "Available columns are:\n",
        paste(nms, collapse = ", ")
      )
    )
  }
  
  return(nms[hit[1]])
}

#########################################################################

### Read GWAS data ###
### sep = "" means whitespace-separated file.
### comment.char = "" is important if the first column is #CHROM.

read_gwas <- function(path, dataset_name) {
  
  df <- read.table(
    file = path,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    fill = TRUE
  )
  
  message("Reading file: ", path)
  message("Number of rows: ", nrow(df))
  message("Number of columns: ", ncol(df))
  message("Columns detected: ", paste(colnames(df), collapse = ", "))
  
  chr_col <- find_col(
    df,
    candidates = c("CHR", "CHROM", "chrom", "chromosome", "#CHROM", "#CHR"),
    required_name = "chromosome"
  )
  
  bp_col <- find_col(
    df,
    candidates = c("BP", "POS", "Position", "position", "GENPOS", "base_pair_location"),
    required_name = "base pair position"
  )
  
  p_col <- find_col(
    df,
    candidates = c("P", "PVAL", "PVALUE", "P_VALUE", "P.value", "p", "pval", "pvalue"),
    required_name = "P value"
  )
  
  df <- df %>%
    rename(
      CHR = all_of(chr_col),
      BP = all_of(bp_col),
      P = all_of(p_col)
    )
  
  df$CHR <- as.character(df$CHR)
  df$BP <- as.numeric(df$BP)
  df$P <- as.character(df$P)
  
  df <- df %>%
    filter(
      !is.na(CHR),
      !is.na(BP),
      !is.na(P),
      P != "",
      P != "NA"
    )
  
  ### Use mpfr only to calculate -log10(P)
  P_mpfr <- mpfr(df$P, precBits = 256)
  
  ### Convert back to normal numeric for dplyr and ggplot
  df$mlog10P <- as.numeric(-log10(P_mpfr))
  
  df <- df %>%
    filter(
      !is.na(mlog10P),
      is.finite(mlog10P)
    )
  
  df$Dataset <- dataset_name
  
  return(df)
}

#########################################################################

### Load Data ###

GWAS_df <- read_gwas(infile, "GWAS_df")
GWAS_df2 <- read_gwas(infile2, "GWAS_df2")

combined_df <- bind_rows(GWAS_df, GWAS_df2)

#########################################################################

### Prepare chromosome order ###

combined_df <- combined_df %>%
  mutate(
    CHR_clean = gsub("^chr", "", CHR, ignore.case = TRUE),
    CHR_num = case_when(
      CHR_clean %in% c("X", "x", "23") ~ 23,
      CHR_clean %in% c("Y", "y", "24") ~ 24,
      CHR_clean %in% c("XY", "xy", "25") ~ 25,
      CHR_clean %in% c("MT", "Mt", "mt", "M", "m", "26") ~ 26,
      TRUE ~ suppressWarnings(as.numeric(CHR_clean))
    )
  ) %>%
  filter(
    !is.na(CHR_num),
    !is.na(BP),
    !is.na(mlog10P)
  ) %>%
  arrange(CHR_num, BP)

#########################################################################

### Create cumulative BP position for x axis ###

cum_BP <- combined_df %>%
  group_by(CHR_clean, CHR_num) %>%
  summarise(
    max_bp = max(BP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CHR_num) %>%
  mutate(
    bp_add = dplyr::lag(cumsum(max_bp), default = 0)
  ) %>%
  select(CHR_clean, CHR_num, bp_add)

Man_plot_df <- combined_df %>%
  inner_join(cum_BP, by = c("CHR_clean", "CHR_num")) %>%
  mutate(
    bp_cum = BP + bp_add,
    plot_y = ifelse(Dataset == "GWAS_df", mlog10P, -mlog10P)
  )

#########################################################################

### Find chromosome centers for x axis labels ###

axis_set <- Man_plot_df %>%
  group_by(CHR_clean, CHR_num) %>%
  summarise(
    center = mean(bp_cum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CHR_num)

#########################################################################

### Set y-axis limit ###

ylim <- Man_plot_df %>%
  filter(Dataset == "GWAS_df") %>%
  summarise(
    ylim = floor(max(mlog10P, na.rm = TRUE)) + 2
  ) %>%
  pull(ylim)

ylim2 <- Man_plot_df %>%
  filter(Dataset == "GWAS_df2") %>%
  summarise(
    ylim = floor(max(mlog10P, na.rm = TRUE)) + 2
  ) %>%
  pull(ylim)

chosen_ylim <- max(ylim, ylim2)

#########################################################################

### Thresholds ###

sig <- p_threshold
sig2 <- p_threshold2

#########################################################################

### Facet labels ###

facet_labels <- c(
  GWAS_df = facetlabel1,
  GWAS_df2 = facetlabel2
)

Man_plot_df$Dataset <- factor(
  Man_plot_df$Dataset,
  levels = c("GWAS_df", "GWAS_df2")
)

#########################################################################

### Miami Plot ###

Miami_plot <- ggplot(
  Man_plot_df,
  aes(
    x = bp_cum,
    y = plot_y,
    color = as_factor(CHR_clean)
  )
) +
  geom_point(size = 0.6, alpha = 0.75) +
  
  geom_hline(
    yintercept = c(-log10(sig), log10(sig)),
    color = "grey40",
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = c(-log10(sig2), log10(sig2)),
    color = "grey60",
    linetype = "dashed"
  ) +
  
  scale_x_continuous(
    labels = axis_set$CHR_clean,
    breaks = axis_set$center
  ) +
  
  scale_color_manual(
    values = rep(
      c("#276FBF", "#183059"),
      length.out = length(unique(axis_set$CHR_clean))
    )
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p-value))
  ) +
  
  facet_wrap(
    ~ Dataset,
    scales = "free_y",
    ncol = 1,
    strip.position = "right",
    labeller = as_labeller(facet_labels)
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 60,
      size = 8,
      vjust = 0.5
    )
  )

#########################################################################

### Set y-axis scales for two facets ###

position_scales <- list(
  
  scale_y_continuous(
    limits = c(0, chosen_ylim),
    breaks = seq(0, chosen_ylim, 1),
    labels = seq(0, chosen_ylim, 1)
  ),
  
  scale_y_continuous(
    limits = c(-chosen_ylim, 0),
    breaks = seq(-chosen_ylim, 0, 1),
    labels = abs(seq(-chosen_ylim, 0, 1))
  )
)

Miami_plot <- Miami_plot +
  facetted_pos_scales(y = position_scales)

#########################################################################

### Save plot ###

outfile <- paste0(filename, ".Miami.pdf")

ggsave(
  plot = Miami_plot,
  filename = outfile,
  width = 40,
  height = 20,
  units = "cm"
)
