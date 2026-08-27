#!/bin/bash

# Download UKB 50k 2.8M sparse shrunk LD matrix for GCTB / SBayesS

### Preamble ###

directory=/media/desk15/iy2120/Project2026/Myproject_AD/02_sex_stratified_gwas/sex_stratified_data/ldm_ukb_50k_bigset_2.8M/

mkdir -p ${directory}
cd ${directory}


### Download ###

wget -c https://zenodo.org/records/3375373/files/ukb_50k_bigset_2.8M.zip.partaa?download=1 -O ukb_50k_bigset_2.8M.zip.partaa

wget -c https://zenodo.org/records/3376357/files/ukb_50k_bigset_2.8M.zip.partab?download=1 -O ukb_50k_bigset_2.8M.zip.partab

wget -c https://zenodo.org/records/3376456/files/ukb_50k_bigset_2.8M.zip.partac?download=1 -O ukb_50k_bigset_2.8M.zip.partac

wget -c https://zenodo.org/records/3376628/files/ukb_50k_bigset_2.8M.zip.partad?download=1 -O ukb_50k_bigset_2.8M.zip.partad

wget -c https://zenodo.org/records/3376628/files/ukb_50k_bigset_2.8M.zip.partae?download=1 -O ukb_50k_bigset_2.8M.zip.partae


### Check md5 ###

md5sum ukb_50k_bigset_2.8M.zip.partaa
md5sum ukb_50k_bigset_2.8M.zip.partab
md5sum ukb_50k_bigset_2.8M.zip.partad
md5sum ukb_50k_bigset_2.8M.zip.partae


### Combine and unzip ###

cat ukb_50k_bigset_2.8M.zip.part* > ukb_50k_bigset_2.8M.zip

unzip ukb_50k_bigset_2.8M.zip