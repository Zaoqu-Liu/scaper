## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(scaper)
library(Seurat)
library(SeuratObject)
library(pheatmap)

## ----message=FALSE, warning=FALSE---------------------------------------------
supportedCytokines(database = "cytosig")
genesetCytoSig(cytokine.eval = "IL6", 
               file.name = system.file("extdata", "IL6_output.csv",
                                       package = "scaper")) %>% head(10)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
# CytoSig.score.output <- scapeForSeurat(seurat.object = pbmc_small,
#                                        database = "cytosig", cytokine = "all", normalize = TRUE)
# class(CytoSig.score.output)
# GetAssay(object = CytoSig.score.output, assay = "scape")
# cytosig_mat <- as.data.frame(t(as.matrix(GetAssayData(CytoSig.score.output, assay = "scape", layer = "data"))))
# pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7,
#          cluster_rows = FALSE, cluster_cols = FALSE)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
# pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
# counts.matrix2 <- as.data.frame(t(as.matrix(GetAssayData(pbmc_small, layer = "data"))))
# CytoSig.score.output <- scape(counts.matrix = counts.matrix2,
#                               database = "cytosig", cytokine = "all")
# pheatmap(CytoSig.score.output, fontsize_row = 4, fontsize_col = 7,
#          cluster_rows = FALSE, cluster_cols = FALSE)

