# Gather cojo_output files, extract gene names, parse each gene output
import sys
import os
import re
import pandas as pd
import parsl
import subprocess
from glob import glob
from subprocess import call
from parsl.app.app import bash_app, python_app, join_app
from parsl.dataflow.futures import AppFuture
from parsl.data_provider.files import File


def parse_gene_gwas(gene_file, gene_name, study, inputs=[], outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0])):
        os.makedirs(os.path.dirname(outputs[0]), exist_ok=True)
    if os.path.exists(outputs[0]):
        return("echo 'Output exists. Remove it or delete it.'")

    output = outputs[0]
    bash_command = \
    f"""
    source activate /media/desk15/iy2120/miniconda3/envs/TWAS
    
    python3 /media/desk15/iy2120/TWAS/Breast-cancer-Example/summary-gwas-imputation/src/gwas_parsing.py \
    -gwas_file {gene_file} \
    -output_column_map SNP variant_id \
    -output_column_map refA effect_allele \
    -output_column_map non_effect_allele non_effect_allele \
    -output_column_map bC effect_size \
    -output_column_map bC_se standard_error \
    -output_column_map Chr chromosome --chromosome_format -output_column_map bp position \
    -output_column_map pC pvalue \
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
    -snp_reference_metadata /media/desk15/iy2120/TWAS/Breast-cancer-Example/data/reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA  \
    -output {output}
    """
    
    proc = subprocess.run(bash_command, shell=True, executable="/bin/bash")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {bash_command}")
    return


cojo_out_files = glob(f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/cojo_output/{study}/*.cma.cojo")
cojo_out_gene_names = [re.match("(ENSG\d{1,}\.\d{1,})\.cma\.cojo", thing)[1] for thing in map(os.path.basename, cojo_out_files)]

gene_parsed_res = []
for gene_file, gene_name in zip(cojo_out_files, cojo_out_gene_names):
	out_file = File(f"/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/condTWAS_parse/{study}/{gene_name}.txt.gz")
	gene_parsed_res.append(parse_gene_gwas(gene_file, gene_name, study=study, inputs=None, outputs=[out_file]))















