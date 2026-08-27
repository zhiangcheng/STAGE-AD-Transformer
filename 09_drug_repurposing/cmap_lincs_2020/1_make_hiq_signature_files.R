library(data.table)
library(dplyr)

############################################################
## Input files
############################################################

## 这是 CMAP / LINCS 的 signature metadata 文件
## 常见名字可能是：
## siginfo_beta.txt
## siginfo_beta.tsv
## siginfo_beta.txt.gz
siginfo_file <- "siginfo_beta.txt"

## 输出目录
setwd("~/Project2026/Myproject_AD/09_drug_repurposing/cmap_lincs_2020")

############################################################
## Read siginfo
############################################################

siginfo <- fread(siginfo_file, sep = "\t", header = TRUE, data.table = FALSE)

cat("Original siginfo dimension:\n")
print(dim(siginfo))

############################################################
## Check required columns
############################################################

required_cols <- c(
  "sig_id",
  "qc_pass",
  "median_recall_rank_spearman",
  "median_recall_rank_wtcs_50"
)

missing_cols <- setdiff(required_cols, colnames(siginfo))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

############################################################
## Optional: only keep compound / drug perturbations
## 如果你的 siginfo 已经对应 level5_beta_trt_cp，
## 这一步不会改变结果。
############################################################

if ("pert_type" %in% colnames(siginfo)) {
  siginfo <- siginfo %>%
    filter(pert_type == "trt_cp")
}

cat("After keeping trt_cp signatures:\n")
print(dim(siginfo))

############################################################
## Clean QC columns
############################################################

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

siginfo$qc_pass_clean <- as.character(siginfo$qc_pass)
siginfo$qc_pass_clean <- trimws(tolower(siginfo$qc_pass_clean))

siginfo$qc_pass_clean <- ifelse(
  siginfo$qc_pass_clean %in% c("1", "true", "t", "yes"),
  1,
  ifelse(
    siginfo$qc_pass_clean %in% c("0", "false", "f", "no"),
    0,
    suppressWarnings(as.numeric(siginfo$qc_pass_clean))
  )
)

siginfo$median_recall_rank_spearman_clean <-
  to_numeric_safe(siginfo$median_recall_rank_spearman)

siginfo$median_recall_rank_wtcs_50_clean <-
  to_numeric_safe(siginfo$median_recall_rank_wtcs_50)

############################################################
## High-quality signature filtering
############################################################

hiq <- siginfo %>%
  filter(
    qc_pass_clean == 1 &
      (
        median_recall_rank_spearman_clean <= 5 |
          median_recall_rank_wtcs_50_clean <= 5
      )
  )

cat("High-quality signatures:\n")
print(dim(hiq))

cat("Number of unique sig_id:\n")
print(length(unique(hiq$sig_id)))

if ("pert_id" %in% colnames(hiq)) {
  cat("Number of unique compounds / perturbagens:\n")
  print(length(unique(hiq$pert_id)))
}

if ("cell_iname" %in% colnames(hiq)) {
  cat("Number of unique cell lines:\n")
  print(length(unique(hiq$cell_iname)))
}

############################################################
## Align high-quality signatures with gctx columns
############################################################

library(cmapR)

gctx_file <- "level5_beta_trt_cp_n720216x12328.gctx"

gctx <- parse.gctx(gctx_file)

gctx_sig_ids <- gctx@cid

hiq <- hiq %>%
  filter(sig_id %in% gctx_sig_ids)

cat("High-quality signatures aligned with gctx:\n")
print(dim(hiq))

cat("Number of unique sig_id after gctx alignment:\n")
print(length(unique(hiq$sig_id)))

############################################################
## Output 1:
## hiq_signature_ids.txt
##
## 这个文件给 3_ensemble_score.R 用：
## ids$sig_id
############################################################

hiq_ids <- data.frame(
  sig_id = unique(hiq$sig_id)
)

write.table(
  hiq_ids,
  file = file.path("hiq_signature_ids.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

############################################################
## Output 2:
## siginfo_beta_hiq_only.txt
##
## 这个文件给 3_ensemble_score.R 用：
## cmap_name, sig_id, pert_id, cell_iname
############################################################

write.table(
  hiq,
  file = file.path("siginfo_beta_hiq_only.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

############################################################
## Optional: also save compressed versions
############################################################

fwrite(
  hiq_ids,
  file = file.path("hiq_signature_ids.txt.gz"),
  sep = "\t",
  quote = FALSE
)

fwrite(
  hiq,
  file = file.path("siginfo_beta_hiq_only.txt.gz"),
  sep = "\t",
  quote = FALSE
)

cat("Done.\n")
cat("Generated files:\n")
cat(file.path("hiq_signature_ids.txt"), "\n")
cat(file.path("siginfo_beta_hiq_only.txt"), "\n")
