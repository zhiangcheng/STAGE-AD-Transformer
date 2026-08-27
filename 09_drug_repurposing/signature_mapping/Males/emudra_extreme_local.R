############################################################
## emudra_extreme_local.R
## Self-contained replacement for:
## eXtremeLogFC.R, eXScores.R, XCor.R, XCorrelation.R,
## XCosine.R, XCos.R, XSum.R
############################################################

## ---------------------------------------------------------
## Helper: clean disease signature
## up_df / down_df should each contain:
##   column 1 = gene symbol
##   column 2 = zscore / logFC / effect direction score
## ---------------------------------------------------------

.prepare_extreme_signature <- function(up_df, down_df) {
  if (is.null(up_df) || nrow(up_df) == 0) {
    up_df <- data.frame(gene = character(), score = numeric())
  }
  
  if (is.null(down_df) || nrow(down_df) == 0) {
    down_df <- data.frame(gene = character(), score = numeric())
  }
  
  up_df <- as.data.frame(up_df)
  down_df <- as.data.frame(down_df)
  
  if (ncol(up_df) < 2) {
    stop("dz_genes_up must have at least 2 columns: gene and score.")
  }
  
  if (ncol(down_df) < 2) {
    stop("dz_genes_down must have at least 2 columns: gene and score.")
  }
  
  up <- data.frame(
    gene = as.character(up_df[[1]]),
    disease_score = as.numeric(as.character(up_df[[2]])),
    direction = "up",
    stringsAsFactors = FALSE
  )
  
  down <- data.frame(
    gene = as.character(down_df[[1]]),
    disease_score = as.numeric(as.character(down_df[[2]])),
    direction = "down",
    stringsAsFactors = FALSE
  )
  
  sig <- rbind(up, down)
  
  sig$gene <- trimws(sig$gene)
  
  sig <- sig[
    !is.na(sig$gene) &
      sig$gene != "" &
      !is.na(sig$disease_score),
  ]
  
  ## If duplicated genes appear, keep the strongest disease signal
  sig <- sig[order(abs(sig$disease_score), decreasing = TRUE), ]
  sig <- sig[!duplicated(sig$gene), ]
  
  rownames(sig) <- NULL
  
  return(sig)
}


## ---------------------------------------------------------
## Generic cosine similarity
## ---------------------------------------------------------

XCosine <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 2) return(NA_real_)
  
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  
  if (denom == 0) return(NA_real_)
  
  sum(x * y) / denom
}


## ---------------------------------------------------------
## Generic correlation
## method = "pearson" or "spearman"
## ---------------------------------------------------------

XCorrelation <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 3) return(NA_real_)
  
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  
  cor(x, y, method = method)
}


## ---------------------------------------------------------
## XSum:
## sum drug scores at disease-up genes
## minus sum drug scores at disease-down genes
##
## Negative XSum usually means drug reverses disease direction:
##   disease-up genes are down in drug
##   disease-down genes are up in drug
## ---------------------------------------------------------

XSum <- function(drug_vec, sig_df) {
  up_genes <- sig_df$gene[sig_df$direction == "up"]
  down_genes <- sig_df$gene[sig_df$direction == "down"]
  
  up_genes <- intersect(up_genes, names(drug_vec))
  down_genes <- intersect(down_genes, names(drug_vec))
  
  up_sum <- if (length(up_genes) > 0) {
    sum(as.numeric(drug_vec[up_genes]), na.rm = TRUE)
  } else {
    0
  }
  
  down_sum <- if (length(down_genes) > 0) {
    sum(as.numeric(drug_vec[down_genes]), na.rm = TRUE)
  } else {
    0
  }
  
  up_sum - down_sum
}


## ---------------------------------------------------------
## XCos:
## cosine similarity between disease scores and drug scores
## on extreme disease genes
## ---------------------------------------------------------

XCos <- function(drug_vec, sig_df) {
  genes <- intersect(sig_df$gene, names(drug_vec))
  
  if (length(genes) < 2) return(NA_real_)
  
  disease_scores <- sig_df$disease_score[match(genes, sig_df$gene)]
  drug_scores <- as.numeric(drug_vec[genes])
  
  XCosine(disease_scores, drug_scores)
}


## ---------------------------------------------------------
## XCor:
## Pearson correlation between disease scores and drug scores
## on extreme disease genes
## ---------------------------------------------------------

XCor <- function(drug_vec, sig_df) {
  genes <- intersect(sig_df$gene, names(drug_vec))
  
  if (length(genes) < 3) return(NA_real_)
  
  disease_scores <- sig_df$disease_score[match(genes, sig_df$gene)]
  drug_scores <- as.numeric(drug_vec[genes])
  
  XCorrelation(disease_scores, drug_scores, method = "pearson")
}


## ---------------------------------------------------------
## XSpe:
## Spearman correlation between disease scores and drug scores
## on extreme disease genes
## ---------------------------------------------------------

XSpe <- function(drug_vec, sig_df) {
  genes <- intersect(sig_df$gene, names(drug_vec))
  
  if (length(genes) < 3) return(NA_real_)
  
  disease_scores <- sig_df$disease_score[match(genes, sig_df$gene)]
  drug_scores <- as.numeric(drug_vec[genes])
  
  XCorrelation(disease_scores, drug_scores, method = "spearman")
}


## ---------------------------------------------------------
## Main wrapper:
##
## mat:
##   rows = genes, rownames(mat) must be gene symbols
##   columns = CMAP signatures / perturbations
##   values = beta / logFC / expression perturbation values
##
## dz_genes_up:
##   col1 = gene symbol
##   col2 = zscore / logFC
##
## dz_genes_down:
##   col1 = gene symbol
##   col2 = zscore / logFC
##
## returns:
##   matrix with columns: xsum, xcos, xcor, xspe
## ---------------------------------------------------------

eXScores <- function(mat, dz_genes_up, dz_genes_down) {
  mat <- as.matrix(mat)
  
  if (is.null(rownames(mat))) {
    stop("mat must have rownames as gene symbols.")
  }
  
  if (is.null(colnames(mat))) {
    stop("mat must have colnames as CMAP signature IDs.")
  }
  
  sig_df <- .prepare_extreme_signature(dz_genes_up, dz_genes_down)
  
  common_genes <- intersect(rownames(mat), sig_df$gene)
  
  if (length(common_genes) < 3) {
    stop(
      paste0(
        "Too few overlapping genes between disease signature and mat: ",
        length(common_genes)
      )
    )
  }
  
  sig_df <- sig_df[sig_df$gene %in% common_genes, ]
  mat_sub <- mat[common_genes, , drop = FALSE]
  
  ## Ensure drug vectors have gene names
  res <- apply(mat_sub, 2, function(drug_vec) {
    names(drug_vec) <- rownames(mat_sub)
    
    c(
      xsum = XSum(drug_vec, sig_df),
      xcos = XCos(drug_vec, sig_df),
      xcor = XCor(drug_vec, sig_df),
      xspe = XSpe(drug_vec, sig_df)
    )
  })
  
  res <- t(res)
  
  return(res)
}