library(data.table)
library(dplyr)
library(biomaRt)

############################################################
## 1. 读取 CMAP row metadata
############################################################

rows <- read.table(
  "./rdesc",
  header = TRUE
)

############################################################
## 2. 提取 Entrez IDs
############################################################

entrez_ids <- unique(as.character(rows$id))

############################################################
## 3. 连接 Ensembl BioMart
############################################################

mart <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  mirror = "uswest"
)

############################################################
## 4. Entrez -> HGNC + Ensembl
############################################################

ids <- getBM(
  attributes = c(
    "entrezgene_id",
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart,
  uniqueRows = TRUE
)

############################################################
## 5. 清理
############################################################

ids <- ids %>%
  filter(
    !is.na(entrezgene_id),
    !is.na(hgnc_symbol),
    !is.na(ensembl_gene_id),
    hgnc_symbol != "",
    ensembl_gene_id != ""
  )

ids$gene <- as.character(ids$entrezgene_id)

ids <- ids %>%
  dplyr::select(
    gene,
    hgnc_symbol,
    ensembl_gene_id
  )

############################################################
## 6. 去重
############################################################

ids <- ids %>%
  distinct(gene, .keep_all = TRUE)

############################################################
## 7. 写出 mat.ids.txt
############################################################

write.table(
  ids,
  file = "mat.ids.txt",
  sep = " ",
  quote = FALSE,
  row.names = FALSE
)

