#' @title Internal VAM Functions
#' @description Internal implementation of the Variance Adjusted Mahalanobis (VAM) method.
#'   These functions are adapted from the VAM package by H. Robert Frost to ensure
#'   compatibility with both Seurat v4 and v5.
#' @author Original: H. Robert Frost <rob.frost@dartmouth.edu>
#' @author Adapted by: Zaoqu Liu <liuzaoqu@163.com>
#' @references Frost (2020) <doi:10.1093/nar/gkaa582>
#' @keywords internal
#' @importFrom MASS fitdistr
#' @importFrom Matrix t
#' @importFrom methods is slotNames

# Core VAM computation function
.vam <- function(gene.expr, tech.var.prop, gene.weights, center = FALSE, gamma = TRUE) {
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }
  
  # Remove genes with 0 mean value
  mean.values <- apply(gene.expr, 2, mean)
  zero.mean.genes <- which(mean.values == 0)
  if (length(zero.mean.genes) > 0) {
    genes.to.keep <- which(mean.values > 0)
    mean.values <- mean.values[genes.to.keep]
    gene.expr <- gene.expr[, genes.to.keep, drop = FALSE]
    if (!missing(tech.var.prop)) {
      tech.var.prop <- tech.var.prop[genes.to.keep]
    }
    if (!missing(gene.weights)) {
      gene.weights <- gene.weights[genes.to.keep]
    }
  }
  
  # Compute gene variances and remove genes with 0 variance
  gene.var <- apply(gene.expr, 2, var)
  zero.var.genes <- which(gene.var == 0)
  if (length(zero.var.genes) > 0) {
    genes.to.keep <- which(gene.var > 0)
    mean.values <- mean.values[genes.to.keep]
    gene.expr <- gene.expr[, genes.to.keep, drop = FALSE]
    gene.var <- gene.var[genes.to.keep]
    if (!missing(tech.var.prop)) {
      tech.var.prop <- tech.var.prop[genes.to.keep]
    }
    if (!missing(gene.weights)) {
      gene.weights <- gene.weights[genes.to.keep]
    }
  }
  
  n <- nrow(gene.expr)
  p <- ncol(gene.expr)
  
  # Use tech.var.prop if specified, otherwise use sample variance
 if (!missing(tech.var.prop)) {
    if (length(tech.var.prop) != p) {
      stop("Length of tech.var.prop not equal to number of genes!")
    }
    tech.var <- tech.var.prop * gene.var
  } else {
    tech.var <- gene.var
  }
  
  # Apply gene weights
  if (!missing(gene.weights)) {
    invalid.weights <- which(gene.weights <= 0)
    if (length(invalid.weights) > 0) {
      stop("gene.weights contains values <= 0!")
    }
    tech.var <- tech.var / gene.weights
  }
  
  # Create inverse covariance matrix
  if (p == 1) {
    inv.cov <- as.matrix(1 / tech.var)
  } else {
    inv.cov <- diag(1 / tech.var)
  }
  
  # Compute squared Mahalanobis distance
  mahalanobis.sq <- apply(gene.expr, 1, function(r) {
    if (center) {
      return((r - mean.values) %*% inv.cov %*% (r - mean.values))
    } else {
      return(r %*% inv.cov %*% r)
    }
  })
  
  if (any(is.infinite(mahalanobis.sq))) {
    stop("Computed infinite distances")
  }
  
  # Compute p-values
  if (gamma) {
    # Gamma distribution method
    gene.expr.null <- apply(gene.expr, 2, function(c) {
      return(sample(c, length(c), replace = FALSE))
    })
    
    mahalanobis.sq.null <- apply(gene.expr.null, 1, function(r) {
      if (center) {
        return((r - mean.values) %*% inv.cov %*% (r - mean.values))
      } else {
        return(r %*% inv.cov %*% r)
      }
    })
    
    nonzero.values <- which(mahalanobis.sq.null > 0)
    gamma.fit <- try(MASS::fitdistr(mahalanobis.sq.null[nonzero.values], "gamma", lower = 0.01), silent = TRUE)
    if (inherits(gamma.fit, "try-error")) {
      p.values <- rep(1, length(mahalanobis.sq))
    } else {
      p.values <- pgamma(mahalanobis.sq, shape = gamma.fit$estimate[1], 
                         rate = gamma.fit$estimate[2], lower.tail = FALSE)
    }
  } else if (center) {
    p.values <- pchisq(mahalanobis.sq, df = p, lower.tail = FALSE)
  } else {
    std.mean.values <- sqrt(mean.values %*% inv.cov %*% mean.values)
    ncp <- sum(std.mean.values^2)
    p.values <- pchisq(mahalanobis.sq, df = p, lower.tail = FALSE, ncp = ncp)
  }
  
  # Ensure no zero p-values
  zero.pvals <- which(p.values == 0)
  if (length(zero.pvals) > 0) {
    min.non.zero <- min(p.values[which(p.values > 0)])
    p.values[zero.pvals] <- min.non.zero
  }
  
  results <- data.frame(cdf.value = 1 - p.values, distance.sq = mahalanobis.sq)
  rownames(results) <- rownames(gene.expr)
  
  return(results)
}

# VAM for multiple gene sets
.vamForCollection <- function(gene.expr, gene.set.collection, tech.var.prop, 
                              gene.weights, center = FALSE, gamma = TRUE) {
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }
  
  p <- ncol(gene.expr)
  n <- nrow(gene.expr)
  
  if (!missing(tech.var.prop)) {
    if (length(tech.var.prop) != p) {
      stop("Length of tech.var.prop does not match number of genes!")
    }
  }
  
  if (!missing(gene.weights)) {
    if (is.list(gene.weights)) {
      if (length(gene.weights) != length(gene.set.collection)) {
        stop("Length of gene.weights list does not match number of gene sets!")
      }
    } else {
      if (length(gene.weights) != p) {
        stop("Length of gene.weights does not match number of genes!")
      }
    }
  } else {
    gene.weights <- rep(1, p)
  }
  
  cell.ids <- rownames(gene.expr)
  num.sets <- length(gene.set.collection)
  set.names <- names(gene.set.collection)
  
  results <- list()
  results$distance.sq <- matrix(0, nrow = n, ncol = num.sets,
                                dimnames = list(cell.ids, set.names))
  results$cdf.value <- matrix(0, nrow = n, ncol = num.sets,
                              dimnames = list(cell.ids, set.names))
  
  for (i in seq_len(num.sets)) {
    set.members <- gene.set.collection[[i]]
    if (length(set.members) == 0) next
    
    set.exprs <- gene.expr[, set.members, drop = FALSE]
    
    if (is.list(gene.weights)) {
      set.weights <- gene.weights[[i]]
    } else {
      set.weights <- gene.weights[set.members]
    }
    
    if (!missing(tech.var.prop)) {
      names(tech.var.prop) <- colnames(gene.expr)
      vam.results <- .vam(gene.expr = set.exprs, tech.var.prop = tech.var.prop[set.members],
                          gene.weights = set.weights, center = center, gamma = gamma)
    } else {
      vam.results <- .vam(gene.expr = set.exprs, gene.weights = set.weights, 
                          center = center, gamma = gamma)
    }
    results$distance.sq[, i] <- vam.results$distance.sq
    results$cdf.value[, i] <- vam.results$cdf.value
  }
  
  return(results)
}

# Helper: Get technical variance proportion for standard normalization (V4)
.getTechVarProp <- function(seurat.data) {
  p <- nrow(seurat.data@assays$RNA@data)
  
  meta.features <- colnames(seurat.data@assays$RNA@meta.features)
  if (length(meta.features) <= 3 ||
      meta.features[2] != "vst.variance" ||
      meta.features[3] != "vst.variance.expected") {
    return(rep(1, p))
  }
  
  tech.var.prop <- seurat.data@assays$RNA@meta.features[, 3] / 
                   seurat.data@assays$RNA@meta.features[, 2]
  tech.var.prop[which(is.nan(tech.var.prop))] <- 1
  
  return(tech.var.prop)
}

# Helper: Get technical variance proportion for V5
.getTechVarPropV5 <- function(seurat.data) {
  # Get number of features
  if ("data" %in% names(seurat.data@assays$RNA@layers)) {
    p <- nrow(seurat.data@assays$RNA@layers$data)
  } else {
    p <- nrow(seurat.data@assays$RNA)
  }
  
  # Check for vst results in meta.data
  if (is.null(seurat.data@assays$RNA@meta.data) || 
      ncol(seurat.data@assays$RNA@meta.data) < 3) {
    return(rep(1, p))
  }
  
  meta.data.cols <- colnames(seurat.data@assays$RNA@meta.data)
  if (!("vf_vst_counts_variance" %in% meta.data.cols) ||
      !("vf_vst_counts_variance.expected" %in% meta.data.cols)) {
    return(rep(1, p))
  }
  
  tech.var.prop <- seurat.data@assays$RNA@meta.data[, "vf_vst_counts_variance.expected"] /
                   seurat.data@assays$RNA@meta.data[, "vf_vst_counts_variance"]
  tech.var.prop[which(is.nan(tech.var.prop))] <- 1
  
  return(tech.var.prop)
}

# Helper: Get technical variance proportion for SCTransform
.getTechVarPropForSCT <- function(seurat.data) {
  SCT.slots <- methods::slotNames(seurat.data@assays$SCT)
  if (!("SCTModel.list" %in% SCT.slots)) {
    stop("Missing slot SCTModel.list!")
  }
  
  num.models <- length(seurat.data@assays$SCT@SCTModel.list)
  if (num.models == 0) {
    stop("SCTModel.list is empty!")
  }
  
  sct.model <- seurat.data@assays$SCT@SCTModel.list[[1]]
  sct.feature.attributes <- names(sct.model@feature.attributes)
  if (!("residual_variance" %in% sct.feature.attributes)) {
    stop("Feature 'residual_variance' missing!")
  }
  
  tech.var.prop <- 1 / sct.model@feature.attributes$residual_variance
  return(tech.var.prop)
}

# Check if Seurat object uses V5 assay
.isAssay5 <- function(assay_obj) {
  inherits(assay_obj, "Assay5")
}

# Get normalized counts from Seurat object (V4/V5 compatible)
.getNormalizedCounts <- function(seurat.data) {
  active_assay <- seurat.data@active.assay
  assay_obj <- seurat.data@assays[[active_assay]]
  
  if (.isAssay5(assay_obj)) {
    # Seurat V5
    if ("data" %in% names(assay_obj@layers)) {
      normalized.counts <- assay_obj@layers$data
      # Ensure row/col names
      if (is.null(rownames(normalized.counts))) {
        rownames(normalized.counts) <- rownames(assay_obj)
      }
      if (is.null(colnames(normalized.counts))) {
        colnames(normalized.counts) <- colnames(seurat.data)
      }
    } else {
      stop("No 'data' layer found in V5 assay!")
    }
  } else {
    # Seurat V4
    normalized.counts <- assay_obj@data
  }
  
  return(normalized.counts)
}

#' Get assay data (V4/V5 compatible)
#' @keywords internal
.getAssayDataCompat <- function(object, assay = NULL, layer = "data") {
  if (is.null(assay)) {
    assay_obj <- object@assays[[object@active.assay]]
  } else {
    assay_obj <- object@assays[[assay]]
  }
  
  if (.isAssay5(assay_obj)) {
    # V5: use layers
    if (layer %in% names(assay_obj@layers)) {
      data_mat <- assay_obj@layers[[layer]]
      if (is.null(rownames(data_mat))) {
        rownames(data_mat) <- rownames(assay_obj)
      }
      if (is.null(colnames(data_mat))) {
        colnames(data_mat) <- colnames(object)
      }
      return(data_mat)
    } else {
      # Try counts as fallback
      return(assay_obj@layers$counts)
    }
  } else {
    # V4: use slots
    if (layer == "data") {
      return(assay_obj@data)
    } else if (layer == "counts") {
      return(assay_obj@counts)
    } else if (layer == "scale.data") {
      return(assay_obj@scale.data)
    }
  }
}

#' VAM for Seurat objects (V4/V5 compatible)
#' @keywords internal
.vamForSeurat <- function(seurat.data, gene.weights, gene.set.collection,
                          center = FALSE, gamma = TRUE, sample.cov = FALSE,
                          return.dist = FALSE) {
  
  if (missing(seurat.data)) {
    stop("Missing Seurat object!")
  }
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }
  
  active_assay <- seurat.data@active.assay
  assay_obj <- seurat.data@assays[[active_assay]]
  is_v5 <- .isAssay5(assay_obj)
  
  if (active_assay == "RNA") {
    normalized.counts <- .getNormalizedCounts(seurat.data)
    
    if (!sample.cov) {
      if (is_v5) {
        tech.var.prop <- .getTechVarPropV5(seurat.data)
      } else {
        tech.var.prop <- .getTechVarProp(seurat.data)
      }
    }
  } else if (active_assay == "SCT") {
    normalized.counts <- .getNormalizedCounts(seurat.data)
    
    if (!sample.cov) {
      tech.var.prop <- .getTechVarPropForSCT(seurat.data)
    }
  } else {
    stop("Unsupported active assay: ", active_assay)
  }
  
  if (missing(gene.weights)) {
    gene.weights <- rep(1, nrow(normalized.counts))
  }
  
  # Compute VAM
  if (sample.cov) {
    vam.results <- .vamForCollection(gene.expr = Matrix::t(normalized.counts),
                                     gene.set.collection = gene.set.collection,
                                     gene.weights = gene.weights,
                                     center = center, gamma = gamma)
  } else {
    vam.results <- .vamForCollection(gene.expr = Matrix::t(normalized.counts),
                                     gene.set.collection = gene.set.collection,
                                     tech.var.prop = tech.var.prop,
                                     gene.weights = gene.weights,
                                     center = center, gamma = gamma)
  }
  
  # Create Assay objects (V4/V5 compatible)
  if (return.dist) {
    dist.mat <- t(vam.results$distance.sq)
    rownames(dist.mat) <- names(gene.set.collection)
    dist.assay <- CreateAssayObject(counts = dist.mat)
    seurat.data[["VAMdist"]] <- dist.assay
  }
  
  cdf.mat <- t(vam.results$cdf.value)
  rownames(cdf.mat) <- names(gene.set.collection)
  cdf.assay <- CreateAssayObject(counts = cdf.mat)
  seurat.data[["VAMcdf"]] <- cdf.assay
  
  return(seurat.data)
}
