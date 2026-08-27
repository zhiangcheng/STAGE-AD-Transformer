rm(list=ls())

tracker <- data.table::fread("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb_tracker.tsv")
gene_list <- gene_list <- stringr::str_extract(tracker$group, "(?<=\\.).+") %>% unique()


# snplist
snplist_file_list <- list.files("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb/snplist/", full.names = TRUE)

library(stringr)
library(data.table)

for (gene in gene_list) {
  # Select files whose names contain the current gene
  matched_files <- snplist_file_list[str_detect(snplist_file_list, gene)]
  
  if (length(matched_files) == 0) {
    warning(paste("No file found for gene:", gene))
    next
  }
  
  # Merge all matching files at once with rbindlist
  df <- rbindlist(lapply(matched_files, fread, header = FALSE))
  df <- unique(df)
  
  # Construct the output path (assuming the same directory and file name = gene.snplist)
  out_path <- file.path("../output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb", 
                        paste0(gene, ".snplist"))
  
  fwrite(df, out_path, col.names = FALSE)
}




# snplist.all.for.prediction
snplist_file_list <- list.files("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb/prediction/", full.names = TRUE)

library(stringr)
library(data.table)

for (gene in gene_list) {
  # Select files whose names contain the current gene
  matched_files <- snplist_file_list[str_detect(snplist_file_list, gene)]
  
  if (length(matched_files) == 0) {
    warning(paste("No file found for gene:", gene))
    next
  }
  
  # Merge all matching files at once with rbindlist
  df <- rbindlist(lapply(matched_files, fread, header = FALSE))
  df <- unique(df)
  
  # Construct the output path (assuming the same directory and file name = gene.snplist)
  out_path <- file.path("../output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb", 
                        paste0(gene, ".snplist.all.for.prediction"))
  
  fwrite(df, out_path, col.names = FALSE)
}




# tracker
head(tracker)

temp <- data.frame(group=gene_list)

temp <- temp %>% mutate(snp_list_path=paste0("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb/",
       group, ".snplist"))


data.table::fwrite(temp,"../output/intermediate_data/study_snp_lists/final_metal_ad_kunkle_pgcalz_ukb_tracker.tsv", 
                   row.names=F, sep="\t")




