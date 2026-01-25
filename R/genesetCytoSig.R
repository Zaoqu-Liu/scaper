#' CytoSig gene set construction.
#' 
#' Returns the gene set associated with input cytokine(s) from the CytoSig database
#' given manually specified cytokine-specific output csv file(s) under 
#' the extdata directory with file name beginning with 
#' the specified cytokine (i.e., 'IL6_output.csv') as currently provided for the IL6 cytokine. 
#' 
#' @param cytokine.eval Cytokine(s) associated with the query.
#' @param file.name List of XML file(s) associated with the cytokine(s) beginning with the specific cytokine name.
#' @return Dataframe consisting of genes and the associated log2 fold change values associated with the specific cytokine(s). 
#' @examples
#' file.name.cytosig1 <- system.file("extdata", "IL6_output.csv", package = "scaper")
#' genesetCytoSig(cytokine.eval = "IL6", file.name = file.name.cytosig1) %>% head(10)
#' @export
genesetCytoSig <- function(cytokine.eval, file.name) {
  # Initialize output container
  output_list <- vector("list", length(cytokine.eval))
  
  # Extract base filenames
  file.name.base <- basename(file.name)
  
  for (j in seq_along(cytokine.eval)) {
    # Match cytokine-specific CSV files
    pattern <- paste0("^", cytokine.eval[j], ".*\\.csv$")
    path.cytokine <- file.name.base[grepl(pattern, file.name.base)]
    
    if (length(path.cytokine) == 0) {
      stop("Specify cytokine-specific output csv file from the CytoSig database with file name beginning with the specified cytokine.")
    }
    
    # Read first matching file
    matched_file <- file.name[grepl(path.cytokine[1], file.name)][1]
    path.out <- read.csv(matched_file, stringsAsFactors = FALSE)
    
    # Ensure numeric conversion
    path.out$Value <- as.numeric(path.out$Value)
    
    # Aggregate log2 fold changes by gene
    path.out2 <- path.out %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(meanVal = mean(Value), .groups = 'drop')
    
    # Construct gene set with weights
    total.out <- data.frame(
      gene = as.character(path.out2$Gene),
      weight = as.numeric(path.out2$meanVal),
      cytokine = cytokine.eval[j],
      stringsAsFactors = FALSE
    )
    
    output_list[[j]] <- total.out
  }
  
  # Combine results across all cytokines
  total.output <- do.call(rbind, output_list)
  
  # Sort by cytokine identifier
  total.output <- total.output %>% dplyr::arrange(cytokine)
  
  return(total.output)
}
