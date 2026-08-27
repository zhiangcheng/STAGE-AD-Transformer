Sys.setLanguage("en")

rm(list=ls())

setwd("/media/desk15/iy2120/Project2026/Myproject_AD/09_drug_repurposing/cmap_lincs_2020/")

library(cmapR)
library(tidyverse)
library(data.table)

gctx <- parse.gctx("level5_beta_trt_cp_n720216x12328.gctx")
sig_info <- fread("siginfo_beta.txt")

mat <- as.data.frame(gctx@mat)
mat <- cbind(gene = gctx@rid, mat)

sig_qc_filtered <- sig_info %>%
  filter(
    qc_pass == 1 &
      (
        median_recall_rank_spearman <= 5 &
          median_recall_rank_wtcs_50 <= 5
      )
  )
valid_sigs <- sig_qc_filtered$sig_id

common_sigs <- intersect(colnames(mat), valid_sigs)
length(duplicated(common_sigs))

mat_qc <- mat[, common_sigs]
mat_qc <- cbind(gene = rownames(mat_qc), mat_qc)

fwrite(mat_qc, "level5_beta_trt_cp_n720216x12328_mat.gz", compress = "gzip")

rdesc <- data.frame(id = gctx@rdesc)
#cdesc <- data.frame(id = gctx@cdesc)
cdesc <- data.frame(id = colnames(mat_qc))

fwrite(rdesc, "rdesc")












