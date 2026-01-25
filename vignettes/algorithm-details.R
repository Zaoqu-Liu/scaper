## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----algorithm, eval=FALSE----------------------------------------------------
# # Pseudocode for bidirectional scoring
# 
# # Step 1: Separate genes by weight sign
# positive_genes <- genes[weights > 0]
# negative_genes <- genes[weights < 0]
# 
# # Step 2: Compute VAM scores
# CDF_positive <- vamForCollection(data, gene.set = positive_genes,
#                                  gene.weights = abs(weights[weights > 0]))
# CDF_negative <- vamForCollection(data, gene.set = negative_genes,
#                                  gene.weights = abs(weights[weights < 0]))
# 
# # Step 3: Bidirectional integration
# Activity <- (CDF_positive + (1 - CDF_negative)) / 2

