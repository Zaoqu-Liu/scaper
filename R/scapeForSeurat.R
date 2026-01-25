#' Cytokine activity scores for a Seurat matrix.
#' 
#' Computes cell-level estimates of cytokine activity for a scRNA-seq Seurat
#' count matrix using the scapeForSeurat method. SCAPE activity estimates are computed by scoring weighted gene 
#' sets from the CytoSig or Reactome databases using the Variance-adjusted Mahalanobis (VAM) method.
#' This function is compatible with both Seurat v4 and v5. Individual gene sets for subsequent scoring 
#' can be reconstructed using the \code{\link[scaper]{genesetCytoSig}} and the \code{\link[scaper]{genesetReactome}} 
#' functions for the CytoSig and the Reactome database, respectively. 
#' @seealso{\code{\link{genesetCytoSig}}, \code{\link{genesetReactome}}, \code{\link{scape}}}
#' 
#' @param seurat.object Seurat counts matrix. 
#' @param database Database used for gene set construction and set scoring.
#' * "cytosig" (default) performs scoring for up to 41 cytokines using the CytoSig database. 
#' * "reactome" performs scoring for up to 30 cytokines using the Reactome database. 
#' @param cytokine Vector of cytokine names to score for activity. The default value of "all" 
#' will score all 41 cytokines supported by CytoSig or 31 supported by Reactome. Please see 
#' function \code{\link[scaper]{supportedCytokines}} to view all the CytoSig or the Reactome 
#' specific scored cytokines. 
#' @param normalize Boolean indicator for whether normalization should be performed before performing gene set scoring. 
#' @return Seurat object consisting of cell-level cytokine activity scores returned as a separate assay (scape for scoring via the CytoSig database and VAMcdf for scoring via the Reactome database). 
#' @examples
#' \donttest{
#' library(SeuratObject)
#' CytoSig.score.output.all <- scapeForSeurat(seurat.object = pbmc_small, 
#' database = "cytosig", cytokine = "all", normalize = TRUE)
#' head(as.data.frame(GetAssayData(CytoSig.score.output.all, assay = "scape", layer = "data")))[1:6,1:3]
#' CytoSig.score.output.specific <- scapeForSeurat(seurat.object = pbmc_small, 
#' database = "cytosig", cytokine = c("IL4", "IL13"), normalize = TRUE)
#' head(as.data.frame(GetAssayData(CytoSig.score.output.specific, assay = "scape", layer = "data")))[,1:3]
#' }
#'
#' @export
scapeForSeurat <- function(seurat.object, database="cytosig", cytokine="all", normalize=TRUE) { 
  if (missing(seurat.object)) {
    stop("scRNA-seq Seurat matrix is missing.")
  }
  
  # Validate database selection
  if (!database %in% c("cytosig", "reactome")) {
    stop("Database must be either 'cytosig' or 'reactome'.")
  }
  
  if (normalize == TRUE) {
    seurat.object <- NormalizeData(seurat.object) 
  }
  
  # Extract gene identifiers from Seurat object
  gene_names <- rownames(seurat.object)
  
  # Validate Seurat object contains genes
  if (length(gene_names) == 0) {
    stop("Seurat object must contain at least one gene.")
  }
  if (database == "cytosig") {
    if (length(cytokine) == 1 && cytokine == "all") {
      cytokine.list <- supportedCytokines(database = "cytosig")
    } else {
      cytokine.list <- cytokine
      # Validate cytokine names
      supported <- supportedCytokines(database = "cytosig")
      invalid <- cytokine.list[!cytokine.list %in% supported]
      if (length(invalid) > 0) {
        stop(paste("Invalid cytokine name(s):", paste(invalid, collapse = ", "),
                   "\nUse supportedCytokines('cytosig') to see valid options."))
      }
    }
    
    # Retrieve cytokine-specific gene sets and weights
    total.output <- mat.cytosig.cytokine.weights[cytokine.list]
    
    # Initialize gene set and weight lists
    list_genes_positive <- vector("list", length(total.output))
    list_weights_positive <- vector("list", length(total.output))
    list_genes_negative <- vector("list", length(total.output))
    list_weights_negative <- vector("list", length(total.output))
    
    names(list_genes_positive) <- paste(cytokine.list, "GenesPositive", sep = "-")
    names(list_weights_positive) <- paste(cytokine.list, "WeightsPositive", sep = "-")
    names(list_genes_negative) <- paste(cytokine.list, "GenesNegative", sep = "-")
    names(list_weights_negative) <- paste(cytokine.list, "WeightsNegative", sep = "-")
    
    for (i in seq_along(total.output)) {
      out.data <- total.output[[i]]
      
      # Separate genes by weight direction
      is_positive <- out.data$weight > 0
      
      # Extract positively weighted genes
      if (any(is_positive)) {
        pos_data <- out.data[is_positive, ]
        gene_match <- pos_data$gene %in% gene_names
        if (any(gene_match)) {
          pos_genes <- pos_data$gene[gene_match]
          list_genes_positive[[i]] <- match(pos_genes, gene_names)
          list_weights_positive[[i]] <- pos_data$weight[gene_match]
        }
      }
      
      # Extract negatively weighted genes
      if (any(!is_positive)) {
        neg_data <- out.data[!is_positive, ]
        gene_match <- neg_data$gene %in% gene_names
        if (any(gene_match)) {
          neg_genes <- neg_data$gene[gene_match]
          list_genes_negative[[i]] <- match(neg_genes, gene_names)
          list_weights_negative[[i]] <- abs(neg_data$weight[gene_match])
        }
      }
    }
    
    # Compute VAM scores for positively weighted genes
    vam.out.positive <- .vamForSeurat(seurat.data = seurat.object, 
                                      gene.weights = list_weights_positive,
                                      gene.set.collection = list_genes_positive)
    
    # Compute VAM scores for negatively weighted genes
    vam.out.negative <- .vamForSeurat(seurat.data = seurat.object, 
                                      gene.weights = list_weights_negative,
                                      gene.set.collection = list_genes_negative)
    
    # Integrate bidirectional gene set scores
    # Negatively weighted genes require complement transformation (1 - CDF)
    # to align directionality with biological interpretation
    pos_data <- .getAssayDataCompat(vam.out.positive, assay = "VAMcdf")
    neg_data <- .getAssayDataCompat(vam.out.negative, assay = "VAMcdf")
    scape.out <- (pos_data + (1 - neg_data))/2
    
    vam.out.positive[['VAMcdf']] <- NULL
    
    # Set row names before creating assay to ensure compatibility
    rownames(scape.out) <- cytokine.list
    scape.assay <- CreateAssayObject(scape.out)
    vam.out.positive[['scape']] <- scape.assay
    
    return(vam.out.positive)
  }
  else {
    if (length(cytokine) == 1 && cytokine == "all") {
      cytokine.list <- supportedCytokines(database = "reactome")
    } else {
      cytokine.list <- cytokine
      # Validate cytokine names
      supported <- supportedCytokines(database = "reactome")
      invalid <- cytokine.list[!cytokine.list %in% supported]
      if (length(invalid) > 0) {
        stop(paste("Invalid cytokine name(s):", paste(invalid, collapse = ", "),
                   "\nUse supportedCytokines('reactome') to see valid options."))
      }
    }
    
    # Retrieve cytokine-specific gene sets
    total.output <- mat.reactome.cytokine[cytokine.list]
    
    # Map gene sets to row indices
    list_genes <- lapply(total.output, function(out.data) {
      gene_match <- out.data$gene %in% gene_names
      if (any(gene_match)) {
        # Use match to preserve gene order
        matched_genes <- out.data$gene[gene_match]
        match(matched_genes, gene_names)
      } else {
        integer(0)
      }
    })
    
    names(list_genes) <- paste(cytokine.list, "Genes", sep = "-")
    
    # Filter out empty gene sets
    list_genes2 <- list_genes[lengths(list_genes) > 0]
    
    if (length(list_genes2) == 0) {
      stop("No genes from the selected cytokines were found in the Seurat object.")
    }
    
    # Compute VAM scores
    vam.out <- .vamForSeurat(seurat.data = seurat.object, 
                             gene.set.collection = list_genes2)
    
    # Update row names using compatible method
    vam_data <- .getAssayDataCompat(vam.out, assay = "VAMcdf")
    
    rownames(vam_data) <- gsub("-Genes", "", rownames(vam_data))
    
    # Recreate assay with corrected names
    vam.out[['VAMcdf']] <- NULL
    vam.out[['VAMcdf']] <- CreateAssayObject(vam_data)
    
    return(vam.out)
  } 
}