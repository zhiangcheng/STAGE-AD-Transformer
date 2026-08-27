library(tidyverse)

STUDY_TAGS <- c("final_metal_ad_kunkle_pgcalz_ukb")
ACAT_TYPES <- c("eqtl", "sqtl")


# TPATH <- glue("../../02_acat_eqtl_sqtl/output/acat_{ACAT_TYPES[1]}/")
# TPATT <- glue("*{STUDY_TAGS[1]}_{ACAT_TYPES[1]}*")

main <- function(study_tag, acat_type, force_overwrite=T){
  require(glue)
  require(tidyverse)
  require(data.table)

  write_name <- glue("../input/{study_tag}_{acat_type}_acat_sig.tsv")

	# `want_files` will take all named tissue sets, the union of distinct genes
	# is run through cojo and then recombined through the tissue list later.
  want_files <- list.files(path=glue("/media/desk15/iy2120/Project2026/Myproject_AD/04_acat_eqtl_sqtl/output/acat_{acat_type}/"),
                           pattern=glue("*{study_tag}_{acat_type}*"),
                           full.names=T)
  print(want_files)

  pval_col <- switch(as.character(acat_type),
                     eqtl="acat",
                     sqtl="mtiss_acat")
  get_sig_gene_col <- function(fpath){
    temp_df <- fread(fpath)
    bfr_thresh <- .05 / temp_df %>%
      filter(!is.na(get(pval_col))) %>%
      nrow()

    rdf <- temp_df %>%
      filter(get(pval_col) <= bfr_thresh) %>%
      select(all_of(c("group")))
      
    return(rdf)
  }
  df_no_dup <- want_files %>%
    map_dfr(~get_sig_gene_col(.)) %>%
    distinct()

  print("before write")
  print(write_name)

  # write_name <- glue("../input/{study_tag}_{acat_type}_acat_sig.tsv")

	if (acat_type == "sqtl"){
		# For sqtl, add brain hippocampus significant genes
		acat_sqtl_breakdown <- fread(glue("/media/desk15/iy2120/Project2026/Myproject_AD/04_acat_eqtl_sqtl/output/acat_tiss_breakdown_sqtl/14tiss__{study_tag}_sqtl_acat_results.tsv")) %>%
			select(group, tissue, tissue_acat) %>%
			filter(tissue == "Brain_Hippocampus") %>%
			distinct() %>%
			select(`group`, `Brain_Hippocampus_sqtl_tissue_acat` = tissue_acat) %>%
			distinct()

		brain_hippocampus_thresh <- .05 / nrow(acat_sqtl_breakdown %>% filter(!is.na(Brain_Hippocampus_sqtl_tissue_acat)))

		acat_sqtl_breakdown <- acat_sqtl_breakdown %>%
			filter(Brain_Hippocampus_sqtl_tissue_acat < brain_hippocampus_thresh)

		print(glue("{study_tag} {acat_type} has  {nrow(df_no_dup)} genes before adding Hippocampus significant only"))
		df_no_dup <- bind_rows(df_no_dup, acat_sqtl_breakdown %>% select(`group`)) %>%
			distinct()

		print(glue("{study_tag} {acat_type} has  {nrow(df_no_dup)} genes after adding Hippocampus significant only"))
	} else if (acat_type == "eqtl"){
		spredixcan_eqtl_hippo <- fread(glue("/media/desk15/iy2120/Project2026/Myproject_AD/03_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_{study_tag}__PM__Brain_Hippocampus.csv"))

		brain_hippo_thresh <- .05 / nrow(spredixcan_eqtl_hippo %>% filter(!is.na(pvalue)))
		spredixcan_eqtl_hippo <- spredixcan_eqtl_hippo %>%
			filter(pvalue < brain_hippo_thresh)
		print(glue("{study_tag} {acat_type} has  {nrow(df_no_dup)} genes before adding Hippocampus significant only"))
		df_no_dup <- bind_rows(df_no_dup, spredixcan_eqtl_hippo %>% select(`group`=`gene`)) %>%
			distinct()
		print(glue("{study_tag} {acat_type} has  {nrow(df_no_dup)} genes after adding Hippocampus significant only"))
	}

	if (file.exists(write_name)){
		check_df <- fread(write_name)
		if (nrow(check_df) != nrow(df_no_dup)){
			if (isTRUE(force_overwrite)){
				write_tsv(x=df_no_dup,
										file=write_name)
			} else {
				message(glue("Not writing {write_name}, file already exists with a differing number of rows."))
				message("")
			}
		} else {
			write_tsv(x=df_no_dup,
									file=write_name)
		}
	} else {
		write_tsv(x=df_no_dup,
								file=write_name)
	}
  return()
}

if (!interactive()){
  run_params_df <- tibble(expand.grid(study_tag=STUDY_TAGS, acat_type=ACAT_TYPES))
  pwalk(run_params_df, main)
} else {
  run_params_df <- tibble(expand.grid(study_tag=STUDY_TAGS, acat_type=ACAT_TYPES))
  pwalk(run_params_df, main)
}
