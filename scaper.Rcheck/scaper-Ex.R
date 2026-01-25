pkgname <- "scaper"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scaper')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("genesetCytoSig")
### * genesetCytoSig

flush(stderr()); flush(stdout())

### Name: genesetCytoSig
### Title: CytoSig gene set construction.
### Aliases: genesetCytoSig

### ** Examples

file.name.cytosig1 <- system.file("extdata", "IL6_output.csv", package = "scaper")
genesetCytoSig(cytokine.eval = "IL6", file.name = file.name.cytosig1) %>% head(10)



cleanEx()
nameEx("genesetReactome")
### * genesetReactome

flush(stderr()); flush(stdout())

### Name: genesetReactome
### Title: Reactome gene set construction.
### Aliases: genesetReactome

### ** Examples

file.name.reactome1 <- system.file("extdata", "IL6_Interleukin6_signaling.xml", 
package = "scaper")
genesetReactome(cytokine.eval = "IL6", file.name = file.name.reactome1) %>% head(10)



cleanEx()
nameEx("scape")
### * scape

flush(stderr()); flush(stdout())

### Name: scape
### Title: Cytokine activity scores for a normalized matrix.
### Aliases: scape

### ** Examples





cleanEx()
nameEx("scapeForSeurat")
### * scapeForSeurat

flush(stderr()); flush(stdout())

### Name: scapeForSeurat
### Title: Cytokine activity scores for a Seurat matrix.
### Aliases: scapeForSeurat

### ** Examples





cleanEx()
nameEx("supportedCytokines")
### * supportedCytokines

flush(stderr()); flush(stdout())

### Name: supportedCytokines
### Title: Gene set to score for the CytoSig or the Reactome databases.
### Aliases: supportedCytokines

### ** Examples

supportedCytokines(database = "cytosig")
supportedCytokines(database = "reactome")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
