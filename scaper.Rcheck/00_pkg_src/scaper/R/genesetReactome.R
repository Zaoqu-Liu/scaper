#' Reactome gene set construction.
#' 
#' Returns the gene set consisting of genes on the Reactome pathway 
#' hierarchy given input cytokine(s) and specified list of pathway 
#' hierarchy xml file(s) with file name beginning with the specified cytokine 
#' (i.e., 'IL6_Interleukin6_signaling.xml') as currently provided for the IL6 cytokine. 
#' 
#' @param cytokine.eval Cytokine(s) associated with the query.
#' @param file.name List of XML file(s) associated with the cytokine(s) beginning with the specific cytokine name. 
#' @return Dataframe consisting of list of genes on the molecular pathway associated with the specified cytokine(s). 
#' @examples
#' file.name.reactome1 <- system.file("extdata", "IL6_Interleukin6_signaling.xml", 
#' package = "scaper")
#' genesetReactome(cytokine.eval = "IL6", file.name = file.name.reactome1) %>% head(10)
#' @export
genesetReactome <- function(cytokine.eval, file.name) {
  # Initialize output container
  output_list <- vector("list", length(cytokine.eval))
  
  # Extract base filenames
  file.name.base <- basename(file.name)
  
  for (j in seq_along(cytokine.eval)){
    # Match cytokine-specific XML pathway files
    pattern <- paste0("^", cytokine.eval[j], ".*\\.xml$")
    filenames <- file.name.base[grepl(pattern, file.name.base)]
    
    if (length(filenames) == 0) {
      stop("Specify Reactome pathway hierarchy xml file(s) with file name beginning with the specified cytokine.")
    }
    
    # Resolve full file paths
    filenames.update <- sapply(filenames, function(fn) {
      file.name[grepl(fn, file.name, fixed = TRUE)][1]
    })
    
    # Parse XML pathway documents
    path.out <- lapply(filenames.update, read_xml)
    
    # Extract gene symbols from each pathway
    str.output.final <- lapply(path.out, function(xml_doc) {
      # Identify biochemical reactions
      allPathwayComponents <- xml_find_all(xml_doc, "/rdf:RDF/bp:biochemicalReaction")
      
      if (length(allPathwayComponents) == 0) {
        return(character(0))
      }
      
      # Extract reaction names
      downstreamSteps <- sapply(seq_along(allPathwayComponents), function(r) {
        stringInquiry <- paste0("/rdf:RDF/bp:biochemicalReaction[", r, "]/bp:NAME")
        xml_text(xml_find_all(xml_doc, stringInquiry))
      })
      
      # Parse gene symbols from text
      downstreamSteps2 <- unlist(downstreamSteps)
      out <- strsplit(downstreamSteps2, "\\s+") 
      out2 <- strsplit(unlist(out), ":")
      str.output <- str_trim(str_extract(unlist(out2), "([[:upper:][0-9]]){3,}"))
      str.output2 <- unique(str.output[!is.na(str.output)])
      
      # Retain alphanumeric gene symbols
      str.output2[grepl("\\D", str.output2)]
    })
    
    # Aggregate genes across pathways
    genes <- unique(unlist(str.output.final))
    
    # Construct gene set
    str.output.final.parsed <- data.frame(
      gene = genes,
      cytokine = cytokine.eval[j],
      stringsAsFactors = FALSE
    )
    
    # Append cognate receptor genes
    receptor.sub <- ramilowski.database %>% 
      dplyr::filter(Ligand.ApprovedSymbol == cytokine.eval[j])
    
    if (nrow(receptor.sub) > 0) {
      receptor.sub1 <- data.frame(
        gene = as.character(receptor.sub$Receptor.ApprovedSymbol),
        cytokine = cytokine.eval[j],
        stringsAsFactors = FALSE
      )
      str.output.final.parsed <- rbind(str.output.final.parsed, receptor.sub1)
    }
    
    output_list[[j]] <- str.output.final.parsed
  }
  
  # Combine results across all cytokines
  str.output.final.parsed.total <- do.call(rbind, output_list)
  
  # Remove duplicate entries and sort
  str.output.final.parsed.total <- unique(str.output.final.parsed.total)
  str.output.final.parsed.total <- str.output.final.parsed.total %>% dplyr::arrange(cytokine)
  
  return(str.output.final.parsed.total)
}
