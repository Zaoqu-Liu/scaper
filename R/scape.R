#' Cytokine activity scores for a normalized matrix.
#' 
#' Computes cell-level estimates of cytokine activity for a normalized scRNA-seq 
#' count matrix using the SCAPE method. SCAPE activity estimates are computed by scoring weighted genes 
#' sets from the CytoSig or Reactome databases using the Variance-adjusted Mahalanobis (VAM) method.
#' Individual gene sets for subsequent scoring can be reconstructed using 
#' the \code{\link[scaper]{genesetCytoSig}} and the \code{\link[scaper]{genesetReactome}} functions 
#' for the CytoSig and the Reactome database, respectively. 
#' @seealso{\code{\link{genesetCytoSig}}, \code{\link{genesetReactome}}}
#' 
#' @param counts.matrix A \eqn{m x n} normalized counts matrix with \eqn{m} samples and \eqn{n} genes. 
#' @param database Database used for gene set construction and set scoring.
#' * "cytosig" performs scoring for up to 41 cytokines using the CytoSig database. 
#' * "reactome" performs scoring for up to 30 cytokines using the Reactome database. 
#' @param cytokine Vector of cytokine names to score for activity. The default value 
#' of "all" will score all 41 cytokines supported by CytoSig or 31 supported by Reactome. Please see 
#' function \code{\link[scaper]{supportedCytokines}} to view all the CytoSig or the Reactome 
#' specific scored cytokines. 
#' @return A \eqn{m x p} matrix consisting of the cell-level cytokine activity scores for \eqn{p} cytokines. 
#' @examples
#' \donttest{
#' library(Seurat)
#' library(SeuratObject)
#' # Note: NormalizeData usage depends on Seurat version
#' pbmc_small <- NormalizeData(pbmc_small)
#' counts.matrix <- as.data.frame(t(as.matrix(GetAssayData(pbmc_small, layer = "data"))))
#' CytoSig.score.output <- scape(counts.matrix = counts.matrix, 
#' database = "cytosig")
#' head(CytoSig.score.output)[,1:3]
#' CytoSig.score.output.specific <- scape(counts.matrix = counts.matrix, 
#' database = "cytosig", cytokine = c("IL4", "IL13"))
#' head(CytoSig.score.output.specific)
#' }
#'
#' @export
scape <- function(counts.matrix, database="cytosig", cytokine="all") { 
  if (missing(counts.matrix)) {
    stop("scRNA-seq counts matrix is missing.")
  }
  
  # Validate database selection
  if (!database %in% c("cytosig", "reactome")) {
    stop("Database must be either 'cytosig' or 'reactome'.")
  }
  
  # Validate matrix dimensions
  if (nrow(counts.matrix) == 0 || ncol(counts.matrix) == 0) {
    stop("Input matrix must contain at least one row and one column.")
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
    
    # Cache column names for gene matching
    col_names <- colnames(counts.matrix)
    
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
        gene_match <- pos_data$gene %in% col_names
        if (any(gene_match)) {
          pos_genes <- pos_data$gene[gene_match]
          list_genes_positive[[i]] <- match(pos_genes, col_names)
          list_weights_positive[[i]] <- pos_data$weight[gene_match]
        }
      }
      
      # Extract negatively weighted genes
      if (any(!is_positive)) {
        neg_data <- out.data[!is_positive, ]
        gene_match <- neg_data$gene %in% col_names
        if (any(gene_match)) {
          neg_genes <- neg_data$gene[gene_match]
          list_genes_negative[[i]] <- match(neg_genes, col_names)
          list_weights_negative[[i]] <- abs(neg_data$weight[gene_match])
        }
      }
    }
    
    # Compute VAM scores for positively weighted genes
    vam.out.positive <- .vamForCollection(gene.expr = counts.matrix, 
                                          gene.weights = list_weights_positive,
                                          gene.set.collection = list_genes_positive)
    
    # Compute VAM scores for negatively weighted genes
    vam.out.negative <- .vamForCollection(gene.expr = counts.matrix, 
                                          gene.weights = list_weights_negative,
                                          gene.set.collection = list_genes_negative)
    
    # Integrate bidirectional gene set scores
    # Negatively weighted genes require complement transformation (1 - CDF)
    # to align directionality with biological interpretation
    vam.out.positive$cdf.value <- (vam.out.positive$cdf.value + (1 - vam.out.negative$cdf.value))/2
    colnames(vam.out.positive$cdf.value) <- cytokine.list
    return(vam.out.positive$cdf.value)
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
    
    # Cache column names for gene matching
    col_names <- colnames(counts.matrix)
    
    # Retrieve cytokine-specific gene sets
    total.output <- mat.reactome.cytokine[cytokine.list]
    
    # Map gene sets to matrix column indices
    list_genes <- lapply(total.output, function(out.data) {
      gene_match <- out.data$gene %in% col_names
      if (any(gene_match)) {
        # Use match to preserve gene order
        matched_genes <- out.data$gene[gene_match]
        match(matched_genes, col_names)
      } else {
        integer(0)
      }
    })
    
    names(list_genes) <- paste(cytokine.list, "Genes", sep = "-")
    
    # Filter out empty gene sets
    list_genes2 <- list_genes[lengths(list_genes) > 0]
    
    if (length(list_genes2) == 0) {
      stop("No genes from the selected cytokines were found in the input matrix.")
    }
    
    # Compute VAM scores
    vam.out <- .vamForCollection(gene.expr = counts.matrix, 
                                 gene.set.collection = list_genes2)
    colnames(vam.out$cdf.value) <- gsub("-Genes", "", 
                                         colnames(vam.out$cdf.value))
    return(vam.out$cdf.value)
  } 
}