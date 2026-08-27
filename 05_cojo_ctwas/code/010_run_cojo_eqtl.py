# -------------------------------COJO------------------------------- #
conda activate GWAS
cd /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/code/
jupyter notebook

# -------------Prepare input data------------- #
import sys
import os
import parsl
import subprocess
import pandas as pd
from subprocess import call
from parsl.app.app import bash_app, python_app, join_app
from parsl.dataflow.futures import AppFuture


AD_INDEX_SNP_WINDOW = 2000000 # 2MB
AD_CLUMP_R2 = 0.65


# Convert the GWAS data to an ma file
def make_cojo_ma_file(study, outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0].filepath)):
        os.makedirs(os.path.dirname(outputs[0].filepath), exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath
    pre_output = os.path.join(os.path.dirname(output), "not_final_" + os.path.basename(output))
    bash_command = \
    f"""
    zcat {GWAS_KEY[study]} | awk 'BEGIN {{FS="\t"}}; {{print $3"_"$4" "$5" "$6" "$8" "$12" "$13" "$11" "$10" {SAMPLE_N_KEY[study]}"}}' > {pre_output}
    """
    proc = subprocess.run(bash_command, shell=True)
    if proc.returncode == 0:
        try:
            pre_cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "zscore", "N"]
            pre_cojo_ma_df = pd.read_csv(pre_output, header=0, names=pre_cojo_ma_names, delim_whitespace=True)

            pre_cojo_ma_df['effect_size_alt'] = pre_cojo_ma_df['zscore'] / np.sqrt(SAMPLE_N_KEY[study])
            pre_cojo_ma_df['N'][pre_cojo_ma_df['N'].isnull()] = SAMPLE_N_KEY[study]

            pre_cojo_ma_df['effect_size'][pre_cojo_ma_df['effect_size'].isnull()] = pre_cojo_ma_df['effect_size_alt'][pre_cojo_ma_df['effect_size'].isnull()]
            pre_cojo_ma_df['standard_error'][pre_cojo_ma_df['standard_error'].isnull()] = 1 / np.sqrt(SAMPLE_N_KEY[study]) 

            cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "N"]
            pre_cojo_ma_df[cojo_ma_names].to_csv(output, sep=" ", header=True, index=False)
            os.remove(pre_output)
        except FileNotFoundError:
            print("Shouldn't have no file.")
    else:
        raise Exception(f"shell command failed? {proc}")
    return

GWAS_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": 
              "/media/desk15/iy2120/Project2026/Myproject_AD/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/final_metal_ad_kunkle_pgcalz_ukb.txt.gz"}
SAMPLE_N_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": 4 / (1/3800 + 1/403806) + 4 / (1/21982 + 1/41944) + 4 / (1/21743 + 1/676035)}
COJO_MA_DIR = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_ma_files/"
study = "final_metal_ad_kunkle_pgcalz_ukb"

cojo_ma_file_res = make_cojo_ma_file(study, outputs=[File(os.path.join(COJO_MA_DIR, f"{study}_gwas.ma"))])

# Create the conditional SNP list file
def make_study_snp_lists(study, gene_whitelist_colname, inputs=[], outputs=[],
                        expression=False,
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    if not os.path.exists(STUDY_SNP_LIST_DIR_KEY[study]):
        os.makedirs(STUDY_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath): 
        return("echo 'Output exists. Remove it or delete it.'")
    cojo_ma_file = inputs[0]
    gene_whitelist  = inputs[1]
    output = outputs[0].filepath

        # python3 /media/desk15/iy2120/TWAS/Breast-cancer-Example/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
    if expression:
        bash_command = \
        f"""
        source activate /media/desk15/iy2120/miniconda3/envs/TWAS
        python3 /media/desk15/iy2120/TWAS/Breast-cancer-Example/summary-gwas-imputation/src/groups_and_conditioned_covariance_for_model.py \
        --gwas_file {GWAS_KEY[study]} \
        --model_db {MODEL_EQTL_DB} \
        --output_dir {STUDY_SNP_LIST_DIR_KEY[study]} \
        --output {output} \
        --parsimony 7 \
        --get_og_for_imputed {cojo_ma_file} \
        --gene_whitelist {gene_whitelist} "{gene_whitelist_colname}"
        """
    else:
        # python3 /media/desk15/iy2120/TWAS/Breast-cancer-Example/MetaXcan/software/SMulTi_gene_snps.py \
        bash_command = \
        f"""
        source activate /media/desk15/iy2120/miniconda3/envs/TWAS
        python3 /media/desk15/iy2120/TWAS/Breast-cancer-Example/MetaXcan/software/SMulTi_gene_snps.py \
        --gwas_file {setup.GWAS_KEY[study]} \
        --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
        --grouping {setup.GROUPING} GTEx_sQTL \
        --model_db_path {setup.MODEL_SQTL_DB} \
        --keep_non_rsid --model_db_snp_key varID \
        --associations {setup.MTISS_TWAS_KEY[study]} SPrediXcan \
        --verbosity 10 \
        --get_og_for_imputed {cojo_ma_file} \
        --snp_list_prefix {setup.STUDY_SNP_LIST_DIR_KEY[study]} \
        --output {output} \
        --gene_whitelist {gene_whitelist} "{gene_whitelist_colname}"
        """

# Execute the command
    result = subprocess.run(bash_command, shell=True, executable='/bin/bash', 
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Command failed. Error details:\n{result.stderr}")
        raise RuntimeError(f"Command failed with code {result.returncode}")
    else:
        print(f"Command completed successfully. Output file: {output}")
        print(result.stdout)
    return result

study = "final_metal_ad_kunkle_pgcalz_ukb"
GENES_TO_RUN = {"final_metal_ad_kunkle_pgcalz_ukb":
						(File("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/final_metal_ad_kunkle_pgcalz_ukb_eqtl_acat_sig.tsv"), "group")}
genes_to_run, colname_for_genes = GENES_TO_RUN[study]
cojo_ma_file_res = File(os.path.join(COJO_MA_DIR, f"{study}_gwas.ma")).path
STUDY_SNP_LIST_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb":
							f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/{study}/"}
STUDY_SNP_LIST_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/"
GWAS_KEY = {"final_metal_ad_kunkle_pgcalz_ukb":
				"/media/desk15/iy2120/Project2026/Myproject_AD/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/final_metal_ad_kunkle_pgcalz_ukb.txt.gz"}
MODEL_EQTL_DB = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/combine_eqtl_for_covar_calc.db"

study_snp_lists_res = make_study_snp_lists(study, gene_whitelist_colname=colname_for_genes, 
										inputs=[cojo_ma_file_res, genes_to_run], outputs=[File(os.path.join(STUDY_SNP_LIST_PARENT, f"{study}_tracker.tsv"))], 
										expression=True)


# Further process the study_snp_lists_res data with prepare_study_snp_lists_eqtl.R


def make_cond_snp_lists(study, inputs=[], outputs=[], bp_range=5000000, imp_check_range=500000, ermatch=None):
                       # stdout=parsl.AUTO_LOGNAME, 
                       # stderr=parsl.AUTO_LOGNAME): 
    def get_gencode(path=GENCODE):
        import pandas as pd
        gencode = pd.read_csv(path, sep="\t", usecols=["chromosome", "start_location", "end_location", "gene_id", "gene_name"])
        gencode["chromosome"] = gencode["chromosome"].str.extract(r"chr(\d{1,2})")
        gencode = gencode.dropna() # Leave out X, Y, M chromosomes (NaN)
        gencode["chr_num"] = gencode["chromosome"].astype(int) # int64 conversion
        return gencode[["chr_num", "start_location", "end_location", "gene_id", "gene_name"]]
    
    def get_snplist(path=REPORTED_SNP_LIST, ermatch=None):
        import numpy as np
        if path == REPORTED_SNP_LIST:
            snplist = pd.read_csv(path, sep="\t")
            snplist = snplist[snplist['p_value'] < 5e-8]
            #if ermatch == "erp":
                #snplist = snplist[snplist['ER+ p'] < 5e-8]
            #elif ermatch == "ern":
                #snplist = snplist[snplist['ER- p'] < 5e-8]
            #else:
                #pass
            snplist['chr_num'] = snplist['chromosome']
            snplist['hg38_SNP'] = snplist.agg('chr{0[chromosome]}_{0[base_pair_location]}'.format, axis=1)
            snplist['Position (hg38)'] = snplist['base_pair_location']

        # hg38_SNP is of the form chr#_<pos>
        return snplist[["chr_num", "hg38_SNP", "Position (hg38)"]]

    import os
    import pandas as pd
    from scipy import stats
    import numpy as np

    if not os.path.exists(COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if not os.path.exists(PROXY_COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(PROXY_COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        print("Output exists. Remove it or delete it.")
        return

    gencode = get_gencode()
    snplist = get_snplist(ermatch=ermatch)

    study_genes = pd.read_csv(inputs[0], sep="\t")
    cojo_ma_df = pd.read_csv(inputs[1], sep=" ", usecols=['chromosome_position', 'effect_size', 'standard_error'])
    cojo_ma_df['position'] = cojo_ma_df['chromosome_position'].str.split("_", expand=True)[1].astype(int)
    cojo_ma_df['chr_num'] = cojo_ma_df["chromosome_position"].str.extract(r"(\d{1,2})").astype(int)
    
    cond_snps_tracker = []
    for gene in study_genes['group']:
        # snplist.all.for.prediction if gene has imputed and original model snps
        # only .snplist if gene only has original model SNPs
        primary_path = "{}".format(os.path.join(STUDY_SNP_LIST_DIR_KEY[study], f"{gene}.snplist.all.for.prediction"
            ))
        secondary_path = "{}".format(os.path.join(STUDY_SNP_LIST_DIR_KEY[study], f"{gene}.snplist"
            ))
        if os.path.exists(primary_path):
            use_path = primary_path
        else:
            use_path = secondary_path

        gene_model_snps = pd.read_csv(use_path, header=None, names=["model_snp"])
        nearby_og_snps_path = 'NA'
        cond_snps_path = 'NA'
        gene_info = gencode[gencode["gene_id"] == gene] # gene_id is ENSG*, gene_name is "WASH7P" sort of stuff
        if len(gene_info) == 0:
            print(f"Can't find this gene in the gencode: {GENCODE}")
            status = "not_found"
            cond_snps_tracker.append((gene, cond_snps_path, None, status, nearby_og_snps_path))
            continue
        elif len(gene_info) > 1:
            print(f"Multiple gencode matches for {gene}, probably repeats, just using the first row.")
            print(gene_info)
            gene_info = gene_info.iloc[[0]]
            # gene_info = gene_info.loc[gene_info.index == 0]
            # status = "multiple_matches" # No index SNPs to condition on
            # cond_snps_tracker.append((gene, cond_snps_path, None, status, nearby_og_snps_path))
            # continue

        gene_chr_num = int(gene_info['chr_num'])
        start_loc = gene_info['start_location'] - bp_range
        end_loc = gene_info['end_location'] + bp_range
        cond_snps = snplist[snplist["Position (hg38)"].between(int(start_loc), int(end_loc))]
        cond_snps = cond_snps[cond_snps['chr_num'] == gene_chr_num]
        cond_snps = cond_snps.drop_duplicates(subset=["hg38_SNP"]) # Remove dups if there for some reason

        if len(cond_snps) == 0:
            print(f"{gene} does't have any reported SNPs within {bp_range} base pairs")
            status = "no_index" # No index SNPs to condition on
            cond_snps_tracker.append((gene, cond_snps_path, gene_chr_num, status, nearby_og_snps_path))
            continue

        # cond_snps = cond_snps.merge(cojo_ma_df[['chromosome_position', 'position', 'effect_size', 'standard_error']], left_on='hg38_SNP', right_on='chromosome_position', how='inner')
        # cond_snps_imp = cond_snps[cond_snps['effect_size'].isna()] # Imputed SNPS
        # cond_snps_og = cond_snps[cond_snps['effect_size'].notna()] # Typed, original SNPs
        cond_snps = cond_snps.merge(cojo_ma_df[["chromosome_position", "position", "effect_size", "standard_error"]], left_on="hg38_SNP", right_on="chromosome_position", how="left")
        cond_snps["position"] = cond_snps["position"].fillna(cond_snps["Position (hg38)"]).astype(int)
        cond_snps_imp = cond_snps[cond_snps["effect_size"].isna()].copy()
        cond_snps_og = cond_snps[cond_snps["effect_size"].notna()].copy()

        # If any of our model snps for the gene overlap with cond_snps_og
        # (original snps used as initial list of conditional snps), save these
        # overlapping SNPs. We will want to 0 them out later after we get COJO
        # results
        gene_model_cond_snps_olap = gene_model_snps.merge(cond_snps_og, left_on='model_snp', 
                right_on='hg38_SNP', how='inner')
        if gene_model_cond_snps_olap.shape[0] > 0:
            model_cond_olap_path = f"{COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist.model.cond.overlap"
            gene_model_cond_snps_olap[["model_snp"]].to_csv(model_cond_olap_path, sep=" ", header=True, index=False)

        if cond_snps_imp.shape[0] > 0:
            # For inspection later, save the original typed SNPs if we're going to add in proxies later for imputed
            cond_snps_path = f"{COND_SNP_LIST_DIR_KEY[study]}{gene}.og.no.proxy.snplist"
            cond_snps_og['zscore'] = cond_snps_og.effect_size / cond_snps_og.standard_error
            cond_snps_og['P'] = 2 * stats.norm.sf(np.abs(cond_snps_og.zscore))
            cond_snps_og = cond_snps_og.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            cond_snps_og[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)

            # Make a list of nearby OG snps
            all_nearby_og_snps = None
            cojo_ma_df_slice = cojo_ma_df[cojo_ma_df['chr_num'] == gene_chr_num]
            cojo_ma_df_slice = cojo_ma_df_slice.dropna() # Want original typed genes only
            for imp_snp in cond_snps_imp.itertuples():
                imp_snp_nearby_og_snps = cojo_ma_df_slice[cojo_ma_df_slice['position'].between(imp_snp.position - imp_check_range - 1, imp_snp.position + imp_check_range + 1)]
                if all_nearby_og_snps is None:
                    all_nearby_og_snps = imp_snp_nearby_og_snps
                else:
                    all_nearby_og_snps = all_nearby_og_snps.append(imp_snp_nearby_og_snps)
                    all_nearby_og_snps = all_nearby_og_snps.drop_duplicates()

            # Don't want to use a model snp for the gene as a potential proxy
            # (low probability, but need to be sure)
            all_nearby_og_snps = all_nearby_og_snps[np.logical_not(
                all_nearby_og_snps['chromosome_position'].isin(gene_model_cond_snps_olap['model_snp']))]

            # Write the list of nearby OG snps
            nearby_og_snps_path = f"{PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.nearbyog.snplist"
            all_nearby_og_snps['zscore'] = all_nearby_og_snps.effect_size / all_nearby_og_snps.standard_error
            all_nearby_og_snps['P'] = 2 * stats.norm.sf(np.abs(all_nearby_og_snps.zscore))
            all_nearby_og_snps = all_nearby_og_snps.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            all_nearby_og_snps[["chromosome_position", "P"]].to_csv(nearby_og_snps_path, sep=" ", header=True, index=False)

            # Write the list of imputed SNPs
            cond_snps_imp_path = f"{PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.imputed.snplist"
            cond_snps_imp['P'] = 0.00000001
            cond_snps_imp[["chromosome_position", "P"]].to_csv(cond_snps_imp_path, sep=" ", header=True, index=False)

            bfile = COJO_BFILE_PATTERN.format(chr_num=gene_chr_num)
            bash_command = f"""
            plink \
              --bfile {bfile} \
              --clump {cond_snps_imp_path},{nearby_og_snps_path} \
              --clump-r2 .7 \
              --clump-kb 500 \
              --clump-p1 0.0000001 \
              --clump-p2 0.000001 \
              --memory 16000 \
              --threads 16 \
              --clump-snp-field chromosome_position \
              --clump-field P \
              --clump-index-first \
              --clump-best \
              --out {PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}
            """
            proc = subprocess.run(bash_command, shell=True)
            if proc.returncode == 0:
                try:
                    og_proxy_snps = pd.read_csv(f"{PROXY_COND_SNP_LIST_DIR_KEY[study]}{gene}.clumped.best", delim_whitespace=True)
                    og_proxy_snps = og_proxy_snps[['PSNP']].astype(str)
                    more_cond_snps_og = og_proxy_snps.merge(cojo_ma_df[['chromosome_position',
                        'position', 'effect_size', 'standard_error']], left_on='PSNP',
                        right_on='chromosome_position', how='inner')
                    more_cond_snps_og = more_cond_snps_og[['chromosome_position', 'position', 'effect_size', 'standard_error']]
                    num_og = len(cond_snps_og)
                    cond_snps_og = pd.concat([cond_snps_og, more_cond_snps_og])
                    num_og_after = len(cond_snps_og)
                    cond_snps_og = cond_snps_og.drop_duplicates()
                    has_imp = True
                except FileNotFoundError: # Sometimes find no clump results
                    pass
                    has_imp = False
            else:
                raise Exception(f"shell command failed? {proc}")
        else:
            has_imp = False


        if len(cond_snps_og) == 0:
            print(f"{gene} does't have any original reported SNPs within {bp_range} base pairs")
            gene_start = int(gene_info['start_location'])
            gene_end = int(gene_info['end_location'])
            min_dist_from_start = abs(cond_snps["Position (hg38)"] - gene_start).min()
            min_dist_from_end = abs(cond_snps["Position (hg38)"] - gene_end).min()
            
            status = "index_snps_imputed" # All index SNPs are imputed in original data, can't run COJO
            if (min_dist_from_start > imp_check_range) or (min_dist_from_end > imp_check_range):
                status += "_closest_over_500_kb_from_gene" 
            cond_snps['P'] = 0.00001

            # # Make a list of nearby OG snps
            # all_nearby_og_snps = None
            # cojo_ma_df_slice = cojo_ma_df[cojo_ma_df.chr_num == gene_chr_num]
            # for imp_snp in cond_snps.itertuples():
            #     index_in_ma = cojo_ma_df_slice.index[cojo_ma_df_slice.chromosome_position == imp_snp.hg38_SNP]
            #     imp_snp_nearby_og_snps = cojo_ma_df_slice.loc[int(index_in_ma.values) - og_check_window:int(index_in_ma.values) + og_check_window]
            #     if all_nearby_og_snps is None:
            #         all_nearby_og_snps = imp_snp_nearby_og_snps
            #     else:
            #         all_nearby_og_snps = all_nearby_og_snps.append(imp_snp_nearby_og_snps)
            #         all_nearby_og_snps = all_nearby_og_snps.drop_duplicates()

            # nearby_og_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.nearbyog.snplist"
            # cond_snps_path = f"{setup.COND_SNP_LIST_DIR_KEY[study]}{gene}.imputed.snplist"
            # cond_snps[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)
            # all_nearby_og_snps['zscore'] = all_nearby_og_snps.effect_size / all_nearby_og_snps.standard_error
            # all_nearby_og_snps['P'] = 2 * stats.norm.sf(np.abs(all_nearby_og_snps.zscore))
            # all_nearby_og_snps = all_nearby_og_snps.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            # all_nearby_og_snps[["hg38_SNP", "P"]].to_csv(nearby_og_snps_path, sep=" ", header=True, index=False)
        else:
            status = "has_original_index_snps"
            if has_imp is True:
                status += "_with_proxies"
            cond_snps_og['zscore'] = cond_snps_og.effect_size / cond_snps_og.standard_error
            cond_snps_og['P'] = 2 * stats.norm.sf(np.abs(cond_snps_og.zscore))
            cond_snps_og = cond_snps_og.sort_values(by='P', ascending=False) # Lowest Pvalue on bottom
            cond_snps_path = f"{COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist"
            cond_snps_og['hg38_SNP'] = cond_snps_og['chromosome_position']
            cond_snps_og[["hg38_SNP", "P"]].to_csv(cond_snps_path, sep=" ", header=True, index=False)
        cond_snps_tracker.append((gene, cond_snps_path, gene_chr_num, status, nearby_og_snps_path))

    tracker_df = pd.DataFrame(cond_snps_tracker, columns=["group", "cond_snp_list_path", "chr_num", "status", "nearby_og_snps_path"])
    tracker_df.to_csv(outputs[0].filepath, sep="\t", index=False)
    return


study = "final_metal_ad_kunkle_pgcalz_ukb"
GENCODE = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/gencode_v26_all.txt"
REPORTED_SNP_LIST = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/Bellenguez_2022_GRCh38.tsv"
COND_SNP_LIST_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/"
COND_SNP_LIST_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb":
							f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/{study}/"}
PROXY_COND_SNP_LIST_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/proxy_work/"
PROXY_COND_SNP_LIST_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb":
							f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/proxy_work/{study}/"}
COJO_BFILE_PATTERN = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/backup_COJO_REF_LD/hg38_EUR_chr_{chr_num}_gtex_v8_1000G.phase3.genotypes.final"
study_snp_lists_res = File(os.path.join(STUDY_SNP_LIST_PARENT, f"{study}_tracker.tsv")).path
cojo_ma_file_res = File(os.path.join(COJO_MA_DIR, f"{study}_gwas.ma")).path

cond_snp_lists_res = make_cond_snp_lists(study, bp_range=AD_INDEX_SNP_WINDOW, 
                                         inputs=[study_snp_lists_res, cojo_ma_file_res], 
                                         outputs=[File(os.path.join(COND_SNP_LIST_PARENT, 
                                                                    f"{study}_tracker.tsv"))], ermatch=None)




def make_clump_cond_snp_list(study, chr_num, plink_output, cond_snp_list,
                      clump_snp_field="hg38_SNP", 
                      clump_field="P",
                      clump_kb=10000,
                      clump_r2=.85,
                      clump_p1=1, # Index variants are chosen greedily, this helps avoid "no clumps made" results
                      cojo_ma_path=None,
                      collin_prevent_mode="plink", # How collinearity issues in COJO --cojo-cond is avoided.
                      outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    if not os.path.exists(CLUMP_COND_SNP_LIST_DIR_KEY[study]):
        os.makedirs(CLUMP_COND_SNP_LIST_DIR_KEY[study], exist_ok=True)
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    bfile = COJO_BFILE_PATTERN.format(chr_num=chr_num)

    final_output = outputs[0]
    out_prefix, _ = os.path.splitext(plink_output)
    plink_bash_command = \
    f"""
    plink \
     --bfile {bfile} \
     --clump {cond_snp_list} \
     --clump-r2 {clump_r2} \
     --clump-kb {clump_kb} \
     --clump-p1 {clump_p1} \
     --memory 16000 \
     --threads 1 \
     --clump-snp-field {clump_snp_field} \
     --clump-field {clump_field} \
     --out {out_prefix}

    if test -f "{plink_output}"; then
      tail -n +2 {plink_output} | tr -s ' ' | cut -d ' ' -f 4 | grep 'chr' > {final_output}
    else
      echo "No plink output? {plink_output}" # Use the most significant old conditioned SNP instead I guess
      cat {cond_snp_list} | tail -n +2 | cut -d ' ' -f 1 | tail -1  > {final_output} 
    fi
    """

    #cojo_ma_path = cojo_ma_path
    gcta_cond_snp_list= cond_snp_list.filepath.replace("snplist", "gcta_snplist")
    gcta_out = out_prefix + ".jma.cojo"
    gcta_final_out = out_prefix + ".clumped.trimmed"

    gcta_bash_command = \
    f"""
    cat {cond_snp_list} | cut -d' ' -f 1 | grep -v 'hg38_SNP' > {gcta_cond_snp_list}


    gcta-1.94.1 \
     --bfile {bfile} \
     --cojo-file {cojo_ma_path} \
     --extract {gcta_cond_snp_list} \
     --cojo-slct \
     --out {out_prefix}

    if test -f "{gcta_out}"; then
      cat {gcta_out} | cut -f 2 | grep -v 'SNP' > {gcta_final_out}
    else
      echo "No SNPs selected. Choosing Top 1"
      cat {gcta_cond_snp_list} | head -1 > {gcta_final_out}
    fi
    """
    if collin_prevent_mode == "gcta":
        bash_cmd = gcta_bash_command
    elif collin_prevent_mode == "plink":
        bash_cmd = plink_bash_command

    proc = subprocess.run(bash_cmd, shell=True, executable="/bin/bash")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {bash_cmd}")
    return final_output


study = "final_metal_ad_kunkle_pgcalz_ukb"
CLUMP_COND_SNP_LIST_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/clump_cond_snp_lists/"
CLUMP_COND_SNP_LIST_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": f"{CLUMP_COND_SNP_LIST_PARENT}{study}/"}
COJO_BFILE_PATTERN = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/backup_COJO_REF_LD/hg38_EUR_chr_{chr_num}_gtex_v8_1000G.phase3.genotypes.final"

cond_snp_lists_res = File(os.path.join(COND_SNP_LIST_PARENT, f"{study}_tracker.tsv")).path
cond_snp_list  = pd.read_csv(cond_snp_lists_res, sep="\t")

clump_cond_snp_lists = []
for gene in cond_snp_list.itertuples():
    if gene.status == "has_original_index_snps":
        clump_cond_snp_lists.append(
                    make_clump_cond_snp_list(study, chr_num=gene.chr_num,
                        plink_output =
                        os.path.join(CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped"),
                        clump_r2 = AD_CLUMP_R2, # Getting collinearity probem with default .85 that worked before
                        cond_snp_list=File(gene.cond_snp_list_path),
                        cojo_ma_path=cojo_ma_file_res,
                        collin_prevent_mode="plink",
                        outputs=[File(os.path.join(CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped.trimmed"))])
            )
    elif gene.status == "has_original_index_snps_with_proxies":
        clump_cond_snp_lists.append(
                    make_clump_cond_snp_list(study, chr_num=gene.chr_num,
                        plink_output =
                        os.path.join(CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped"),
                        clump_r2 = AD_CLUMP_R2, # Getting collinearity probem with default .85 that worked before
                        cond_snp_list=File(gene.cond_snp_list_path),
                        cojo_ma_path=cojo_ma_file_res,
                        collin_prevent_mode="plink",
                        outputs=[File(os.path.join(CLUMP_COND_SNP_LIST_DIR_KEY[study],
                            f"{gene.group}.clumped.trimmed"))])
            )
        # elif gene.status.startswith("has_index_snps_imputed"):
        #     clump_cond_snp_lists.append(
        #             make_clump_og_proxy_cond_snp_list(study, chr_num=gene.chr_num,
        #                 plink_output =
        #                 os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
        #                     f"{gene.group}.clumped.best"),
        #                 cond_snp_list=File(gene.cond_snp_list_path),
        #                 nearby_og_snp_list=File(gene.nearby_og_snps_path),
        #                 outputs=[File(os.path.join(setup.CLUMP_COND_SNP_LIST_DIR_KEY[study],
        #                     f"{gene.group}.clumped.trimmed"))])
        #     )
    else:
        clump_cond_snp_lists.append(None)




def make_cojo_run_data(study, study_list, cond_genes_list, gene_col="group", sep="\t", inputs=[], outputs=[]):
    
    if not os.path.exists(COJO_IN_DIR_KEY[study]):
        os.makedirs(COJO_IN_DIR_KEY[study], exist_ok=True)

    clump_cond_snp_lists = inputs
    whitelist = pd.read_csv(study_list, sep=sep) # group
    cond_genes_df = pd.read_csv(cond_genes_list, sep=sep)
    cond_genes_df['clump_cond_snp_list_path'] = clump_cond_snp_lists
    whitelist = whitelist.merge(cond_genes_df, how='left') # cond_snp_list_path chr_num

    blacklist = whitelist[whitelist['clump_cond_snp_list_path'].isna()] # If no conditional SNPs then can't run COJO
    whitelist = whitelist[whitelist['clump_cond_snp_list_path'].notna()] 

    blist_path, wlist_path = outputs 
    blacklist.to_csv(blist_path, index=False, sep="\t")
    whitelist.to_csv(wlist_path, index=False, sep="\t")
    return

COJO_IN_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_input/"
COJO_IN_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": f"{COJO_IN_PARENT}{study}/"}

cojo_run_dat_res = make_cojo_run_data(study, study_list=study_snp_lists_res,
                                          cond_genes_list=cond_snp_lists_res,
                                          inputs = [elmnt.path 
                                                    if elmnt is not None else None 
                                                    for elmnt in clump_cond_snp_lists],
                                          outputs=[File(os.path.join(COJO_IN_PARENT, f"{study}_blacklist.tsv")),
                                                   File(os.path.join(COJO_IN_PARENT, f"{study}_whitelist.tsv"))])




# Kick off cojo runs
cojo_run_dat_res = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_input/final_metal_ad_kunkle_pgcalz_ukb_whitelist.tsv"
gene_whitelist  = pd.read_csv(cojo_run_dat_res, sep="\t")




# Run COJO
def gcta_cojo(study, gene_name, snp_list_path, cond_snp_list_path, chr_num,
              out_prefix, inputs=[], outputs=[], stdout=parsl.AUTO_LOGNAME,
              stderr=parsl.AUTO_LOGNAME): 
    
    if not os.path.exists(COJO_OUT_DIR_KEY[study]):
        os.makedirs(COJO_OUT_DIR_KEY[study], exist_ok=True)
    
    bfile = COJO_BFILE_PATTERN.format(chr_num=chr_num)
    cojo_ma_file = inputs[0]
    output = outputs[0]
    
    extract_path = "{}{}.snplist".format(COJO_IN_DIR_KEY[study], gene_name)
    pre_cojo = f"""
    IN_DIR=`dirname {extract_path}`
    mkdir -p $IN_DIR
    cat {snp_list_path} {cond_snp_list_path} | sort --unique > {extract_path}
    OUT_DIR=`dirname {out_prefix}`
    mkdir -p $OUT_DIR
    """

    bash_command = \
    f"""
    if test -f "{output}"; then
      echo "Output exists: {output}"
    else
      {pre_cojo}

      gcta-1.94.1 \
      --bfile {bfile} \
      --cojo-file {cojo_ma_file} \
      --cojo-cond {cond_snp_list_path} \
      --extract {extract_path} \
      --cojo-slct \
      --out {out_prefix}
    fi
    """

    proc = subprocess.run(bash_command, shell=True, executable="/bin/bash")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {bash_command}")
    return output

def add_non_effect_allele(study, cojo_ma_path, cojo_res_sep="\t", cojo_ma_sep=" ", inputs=[]):
    import numpy as np
    import pandas as pd 
    from glob import glob
    import os
    import re
    cojo_ma_df = pd.read_csv(cojo_ma_path, sep=cojo_ma_sep)

    cojo_oput_path = inputs[0]
    oput_df = pd.read_csv(cojo_oput_path, sep=cojo_res_sep)
    gene = re.split(".cma.cojo", os.path.basename(cojo_oput_path))[0] # Could be a "protein" or "sequence id" maybe soon too

    # Put in non effect allele if not present
    if 'non_effect_allele' not in oput_df:
        oput_df['conditional_status'] = None
        # If there were any SNPs used to model the gene that overlapped with the original index SNP list, join them in
        # if necessary and zero out their effect size/standard error
        model_cond_olap_path = f"{COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist.model.cond.overlap"
        if os.path.exists(model_cond_olap_path):
            model_cond_olap_df = pd.read_csv(model_cond_olap_path)
            model_cond_olap_df.rename(columns={"model_snp": "SNP"}, inplace=True)
            oput_df = oput_df.merge(model_cond_olap_df, left_on='SNP', right_on='SNP', how='outer')

            # Only zero out conditional overlaps if we otherwise failed to calculate the effect
            # E.g. there is overlap, but the overlapped SNP is dropped of the condition set in clumping
            # and we were able to calculate conditional effects. So we don't throw the effect away.
            oput_df['pC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['pC'].isnull())] = 1
            oput_df['bC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC'].isnull())] = 0
            oput_df['conditional_status'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = "null_result_and_condition_overlap"
            oput_df['bC_se'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = 0

        # Add in any predictors in the model that were imputed. 0 them out here as well.
        # With the zscore and standard error estimation there shouldn't be any imputed SNPs remaining
        gene_all_pred_snps_path = f"{STUDY_SNP_LIST_DIR_KEY[study]}{gene}.snplist.all.for.prediction"
        if os.path.exists(gene_all_pred_snps_path):
            model_all_pred_snps_df = pd.read_csv(gene_all_pred_snps_path, names=["SNP"])
            oput_df = oput_df.merge(model_all_pred_snps_df, left_on='SNP', right_on='SNP', how='outer')

        oput_df = oput_df.merge(cojo_ma_df, left_on='SNP', right_on='chromosome_position')

        # If we inserted any new variant, make sure and Chr, refA
        # (effect allele), bp, and freq are added.
        # Get chromosome
        oput_df['Chr'] = oput_df['SNP'][0].split("_")[0].split("chr")[1]
        oput_df['refA'] = oput_df['effect_allele']
        oput_df['freq'] = oput_df['frequency']

        # For any NaN conditional effects (bC, BC_se, pC), set them to 0, 0, 1
        oput_df['pC'][oput_df['pC'].isnull()] = 1
        oput_df['bC'][oput_df['bC'].isnull()] = 0
        oput_df['conditional_status'][oput_df['bC_se'].isnull()] = "freq_mismatch_or_not_in_ld_or_collinear"
        oput_df['bC_se'][oput_df['bC_se'].isnull()] = 0
        oput_df['conditional_status'][oput_df['conditional_status'].isnull()] = "cojo"

        oput_df[["throwaway_chr", "bp"]] = oput_df['SNP'].str.split("_", expand=True)
        oput_df = oput_df[['Chr', 'SNP', 'bp', 'refA', 'non_effect_allele', 'freq', 'b', 'se', 'p', 'n', 'freq_geno', 'bC', 'bC_se', 'pC', 'conditional_status']]

        new_out = f"{cojo_oput_path}"
        oput_df.to_csv(new_out, sep=cojo_res_sep, index=False)
    return

COJO_IN_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_input/"
COJO_IN_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": f"{COJO_IN_PARENT}{study}/"}
COJO_OUT_PARENT = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_output/"
COJO_OUT_DIR_KEY = {"final_metal_ad_kunkle_pgcalz_ukb": f"{COJO_OUT_PARENT}{study}/"}

cojo_all_res = []
for gene in gene_whitelist.itertuples():
        out_prefix = "{}{}".format(COJO_OUT_DIR_KEY[study], gene.group)
        cojo_all_res.append(gcta_cojo(study, gene.group, gene.snp_list_path,
                  gene.clump_cond_snp_list_path, gene.chr_num, out_prefix,
                  inputs=[cojo_ma_file_res],
                  outputs=[File(f"{out_prefix}.cma.cojo")])
                  )

ne_res = []
for cojo_res in cojo_all_res_path:
    if cojo_res is None:
        print("Skip: cojo_res is None")
        continue
    if not os.path.exists(cojo_res):
        print(f"Skip missing file: {cojo_res}")
        continue

    ne_res.append(
        add_non_effect_allele(
            study,
            cojo_ma_path=cojo_ma_file_res,
            inputs=[cojo_res]
        )
    )


# No longer needed
# In some cases we might have {gene}.cma.allna files, indicating that all the genes we wanted conditional effects for could not be calculated (usually when the # of genes we want is small), so we go in and "fix" 
# these genes by creating a new list of SNPs, taking the nearest <na_snp_window> original SNPs on each side of the NA SNPs and running COJO with those instead.
# make_alt_snp_lists(study, na_snp_window, whitelist_path=cojo_run_dat_res.outputs[1].result().filepath,
#         cojo_ma_path=cojo_ma_file_res.outputs[0].result().filepath)
# parsl.wait_for_current_tasks()
# gene_whitelist  = pd.read_csv(cojo_run_dat_res.outputs[1].result(), sep="\t")
# for gene in gene_whitelist.itertuples():
#     if gene.status == "has_original_index_snps_na_fix":
#         out_prefix = "{}{}".format(COJO_OUT_DIR_KEY[study], gene.group)
#         gcta_cojo(study, gene.group, gene.na_fix_snp_list_path,
#                   gene.clump_cond_snp_list_path, gene.chr_num, out_prefix,
#                   inputs=[cojo_ma_file_res.outputs[0]],
#                   outputs=[File(f"{out_prefix}.cma.cojo")])
# parsl.wait_for_current_tasks()
# add_non_effect_allele(study, cojo_ma_path=cojo_ma_file_res.outputs[0])




# Add APOE(ENSG00000130203.9) gene Since Charlie might want it and run COJO for it
def special_gcta_cojo(study, inputs=[], outputs=[], stdout=parsl.AUTO_LOGNAME,
              stderr=parsl.AUTO_LOGNAME): 
    
    bash_command = \
    f"""
    cp /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/{study}/ENSG00000130203.9.snplist /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/clump_cond_snp_lists/{study}/ENSG00000130203.9.clumped.trimmed

    gcta-1.94.1 \
    --bfile /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/backup_COJO_REF_LD/hg38_EUR_chr_19_gtex_v8_1000G.phase3.genotypes.final \
    --cojo-file /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_ma_files/{study}_gwas.ma \
    --cojo-cond /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cond_snp_lists/{study}/ENSG00000130203.9.snplist \
    --extract /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/study_snp_lists/{study}/ENSG00000130203.9.snplist.all.for.prediction \
    --out /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_output/{study}/APOE/ENSG00000130203.9
    """
    
    proc = subprocess.run(bash_command, shell=True, executable="/bin/bash")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {bash_command}")
    return

def special_add_non_effect_allele(study, cojo_res_sep="\t", cojo_ma_sep=" ", inputs=[],
        outputs=[], stdout=parsl.AUTO_LOGNAME, stderr=parsl.AUTO_LOGNAME):
    import numpy as np
    import pandas as pd
    from glob import glob
    import os
    import re
    cojo_ma_df = pd.read_csv("/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_ma_files/" + study + "_gwas.ma", sep=cojo_ma_sep)

    cojo_oput_path = "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_output/" + study + "/ENSG00000130203.9.cma.cojo"
    oput_df = pd.read_csv(cojo_oput_path, sep=cojo_res_sep)
    gene = re.split(".cma.cojo", os.path.basename(cojo_oput_path))[0] # Could be a "protein" or "sequence id" maybe soon too

    # Put in non effect allele if not present
    if 'non_effect_allele' not in oput_df:
        print(study)
        print(gene)
        oput_df['conditional_status'] = None
        # If there were any SNPs used to model the gene that overlapped with the original index SNP list, join them in
        # if necessary and zero out their effect size/standard error
        model_cond_olap_path = f"{COND_SNP_LIST_DIR_KEY[study]}{gene}.snplist.model.cond.overlap"
        if os.path.exists(model_cond_olap_path):
            model_cond_olap_df = pd.read_csv(model_cond_olap_path)
            model_cond_olap_df.rename(columns={"model_snp": "SNP"}, inplace=True)
            oput_df = oput_df.merge(model_cond_olap_df, left_on='SNP', right_on='SNP', how='outer')

            # Only zero out conditional overlaps if we otherwise failed to calculate the effect
            # E.g. there is overlap, but the overlapped SNP is dropped of the condition set in clumping
            # and we were able to calculate conditional effects. So we don't throw the effect away.
            oput_df['pC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['pC'].isnull())] = 1
            oput_df['bC'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC'].isnull())] = 0
            oput_df['conditional_status'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = "null_result_and_condition_overlap"
            oput_df['bC_se'][np.logical_and(oput_df['SNP'].isin(model_cond_olap_df['SNP']), oput_df['bC_se'].isnull())] = 0

        # Add in any predictors in the model that were imputed. 0 them out here as well.
        # With the zscore and standard error estimation there shouldn't be any imputed SNPs remaining
        gene_all_pred_snps_path = f"{STUDY_SNP_LIST_DIR_KEY[study]}{gene}.snplist.all.for.prediction"
        if os.path.exists(gene_all_pred_snps_path):
            model_all_pred_snps_df = pd.read_csv(gene_all_pred_snps_path, names=["SNP"])
            oput_df = oput_df.merge(model_all_pred_snps_df, left_on='SNP', right_on='SNP', how='outer')

        oput_df = oput_df.merge(cojo_ma_df, left_on='SNP', right_on='chromosome_position')

        # If we inserted any new variant, make sure and Chr, refA
        # (effect allele), bp, and freq are added.
        # Get chromosome
        oput_df['Chr'] = oput_df['SNP'][0].split("_")[0].split("chr")[1]
        oput_df['refA'] = oput_df['effect_allele']
        oput_df['freq'] = oput_df['frequency']

        # For any NaN conditional effects (bC, BC_se, pC), set them to 0, 0, 1
        oput_df['pC'][oput_df['pC'].isnull()] = 1
        oput_df['bC'][oput_df['bC'].isnull()] = 0
        oput_df['conditional_status'][oput_df['bC_se'].isnull()] = "freq_mismatch_or_not_in_ld_or_collinear"
        oput_df['bC_se'][oput_df['bC_se'].isnull()] = 0
        oput_df['conditional_status'][oput_df['conditional_status'].isnull()] = "cojo"

        oput_df[["throwaway_chr", "bp"]] = oput_df['SNP'].str.split("_", expand=True)
        oput_df = oput_df[['Chr', 'SNP', 'bp', 'refA', 'non_effect_allele', 'freq', 'b', 'se', 'p', 'n', 'freq_geno', 'bC', 'bC_se', 'pC', 'conditional_status']]

        new_out = f"{cojo_oput_path}"
        written = oput_df.to_csv(new_out, sep=cojo_res_sep, index=False)
        print("written")
        print(" ")
    return(written)

if study in ["final_metal_ad_kunkle_pgcalz_ukb"]:
        # meta_analysis_BCAC_UKB_ovr_brca has the file but still doesn't run, overwrite
        # Check if this is still needed before running it blindly, only need to
        # run if we didn't have anything to condition this gene on
        # if not os.path.exists(f"../output/intermediate_data/clump_cond_snp_lists/{study}/ENSG00000130203.9.clumped.trimmed"):
        special_res = special_gcta_cojo(study, outputs=[File(f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_output/{study}/APOE/ENSG00000130203.9.cma.cojo")])
        special_add_non_effect_allele(study=study, inputs=[special_res])












