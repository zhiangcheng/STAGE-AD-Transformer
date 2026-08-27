
## code for example tissue (GTEx-ACC)

library(tidyverse)
library(data.table)
library(plyr)
library(dplyr)

#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/eXtremeLogFC.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/eXScores.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/XCor.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/XCorrelation.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/XCosine.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/XCos.R")
#source("/sc/arion/projects/psychgen/alanna/sigmap_prelim/emudra/Source/EMUDRA/R/XSum.R")
source("emudra_extreme_local.R")

args <- commandArgs(trailingOnly=TRUE)
number=args[1]

### calculate eXtreme scores only ###
mat <- fread(paste0('/media/desk15/iy2120/Project2026/Myproject_AD/09_drug_repurposing/cmap_lincs_2020/level5_beta_trt_cp_n720216x12328_mat.gz'))
mat <- as.data.frame(mat)
rows <- read.table("/media/desk15/iy2120/Project2026/Myproject_AD/09_drug_repurposing/cmap_lincs_2020/rdesc",header=TRUE)

library(biomaRt)
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",mirror="uswest")
Ensemble2HGNC <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                       filters = "entrezgene_id", values = rows$id,
                       mart = mart, uniqueRows = TRUE)
Ensemble2HGNC <- Ensemble2HGNC[complete.cases(Ensemble2HGNC),]
Ensemble2HGNC <- Ensemble2HGNC[!duplicated(Ensemble2HGNC$hgnc_symbol),]

Ensemble2HGNC$gene <- Ensemble2HGNC$hgnc_symbol
Ensemble2HGNC$entrezgene_id <- as.character(as.numeric(Ensemble2HGNC$entrezgene_id))

rownames(mat) <- rows$id
#mat <- mat %>% rownames_to_column('gene')
mat$gene <- as.character(mat$gene)
mat <- mat %>% left_join(Ensemble2HGNC, by=c("gene" = "entrezgene_id"))
mat <- mat[,-c(1,(ncol(mat)-1))]

mat <- mat %>% drop_na(gene.y)
mat <- mat %>%
  filter(!is.na(gene.y) & gene.y != "")
rownames(mat) <- mat$gene.y
mat$gene.y <- NULL
# 11,349 genes included in cmap resource after converting from entrez to hgnc id


################################################
######## TWAS QVAL #############################
seq <- c('0.01','0.05','0.1','0.2')

for (num in seq){
    twas <- read.table(paste0('/media/desk15/iy2120/Project2026/Myproject_AD/03_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_final_metal_ad_kunkle_pgcalz_ukb__PM__Brain_Cortex.csv'),header=T,sep=',')
    colnames(twas) <- c('gene','gene_name','zscore','effect_size','pvalue','VAR_g','pred_perf_R2','pred_perf_p','pred_perf_q','n_snps_used','n_snps_in_cov','n_snps_in_model','qval','weight')
    twas$up_down[as.numeric(as.character(twas$zscore))>0] <- "up"
    twas$up_down[as.numeric(as.character(twas$zscore))<0] <- "down"
    dz_signature <- twas
    dz_genes_up <- subset(dz_signature,up_down=="up",select=c("gene_name","zscore"))
    dz_genes_down <- subset(dz_signature,up_down=="down",select=c("gene_name","zscore"))

    x.scores.qval.dat <- eXScores(mat, dz_genes_up, dz_genes_down)
    colnames(x.scores.qval.dat) <- c("xsum","xcos","xcor","xspe")
    x.scores.qval.dat <- as.data.frame(x.scores.qval.dat) %>% rownames_to_column('signature')

    qval.long <- gather(x.scores.qval.dat,score_type,value,xsum:xspe)
    qval.long$twas_sig <- num
    assign(paste("qval.long",num,sep=""),qval.long)
}


################################################
######## TWAS FINEMAPPED ########################
twas <- read.table('/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/output/final_output/focus_eqtl_14tiss_final_metal_ad_kunkle_pgcalz_ukb_by_tissue.focus.tsv',sep="\t",header=T)
#twas <- twas[!(twas$block=="block"),]
twas <- twas[!(twas$mol_name == "" | is.na(twas$mol_name)), ]

twas$twas_z_pop1 <- as.numeric(as.character(twas$twas_z))
twas$up_down[twas$twas_z_pop1>0] <- "up"
twas$up_down[twas$twas_z_pop1<0] <- "down"
dz_signature <- twas
dz_genes_up <- subset(dz_signature,up_down=="up",select=c("mol_name","twas_z_pop1"))
dz_genes_down <- subset(dz_signature,up_down=="down",select=c("mol_name","twas_z_pop1"))

x.scores.finemap.dat <- eXScores(mat, dz_genes_up, dz_genes_down)
colnames(x.scores.finemap.dat) <- c("xsum","xcos","xcor","xspe")
x.scores.finemap.dat <- as.data.frame(x.scores.finemap.dat) %>% rownames_to_column('signature')

finemap.long <- gather(x.scores.finemap.dat,score_type,value,xsum:xspe)
finemap.long$twas_sig <- "credset"


master <- bind_rows(qval.long0.01,qval.long0.05,qval.long0.1,qval.long0.2,finemap.long)

write.table(master,paste0("/media/desk15/iy2120/Project2026/Myproject_AD/09_drug_repurposing/signature_mapping/Females/extreme_scores_acc_femalead.txt"),quote=F,row.names=F)

