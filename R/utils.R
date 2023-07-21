## identify column names
identifyColNames <- function(tabOutput, tabID = "") {

  # define accepted column names
  idAccepted <- "ID|Symbol|Gene_Symbol|Mirna|mir|Gene|gene.symbol|Gene.symbol"
  fcAccepted <- "logFC|log2FoldChange|FC|lFC"
  exprAccepted <- "AveExpr|baseMean|logCPM"
  pvalAccepted <- "P.Value|pvalue|PValue|Pvalue"
  fdrAccepted <- "adj.P.Val|padj|FDR|fdr|adj|adj.p|adjp"
  acceptedNames <- c(idAccepted, fcAccepted, exprAccepted,
                     pvalAccepted, fdrAccepted)

  ## try to identify column names
  dfNames <- colnames(tabOutput)
  idCol <- grep(idAccepted, dfNames)
  fcCol <- grep(fcAccepted, dfNames)
  exprCol <- grep(exprAccepted, dfNames)
  pvalCol <- grep(pvalAccepted, dfNames)
  fdrCol <- grep(fdrAccepted, dfNames)
  tableCols <- list(idCol, fcCol, exprCol, pvalCol, fdrCol)
  tableNames <- c("ID", "logFC", "AveExpr", "P.Value", "adj.P.Val")

  ## check if columns are correctly identified
  if (any(lengths(tableCols) == 0)) {
    stop(paste("The software is unable to automatically find columns",
               "relative to",
               paste(tableNames[which(lengths(tableCols) == 0)],
                     collapse = ", "),
               "in the", tabID, "differential expression data.frame!",
               "Please change the name of these columns to one of:",
               paste(acceptedNames[which(lengths(tableCols) == 0)],
                     collapse = ", ")), call. = FALSE)
  } else if (any(lengths(tableCols) > 1)) {
    stop(paste("More than one column can be interpreted as",
               paste(tableNames[which(lengths(tableCols) > 1)],
                     collapse = ", "),
               "in the", tabID, "differential expression data.frame!",
               "Please rename these columns to have unambiguous column names!",
               "See ?MirnaExperiment for further details"), call. = FALSE)
  }

  ## create a data.frame with desired columns
  tabOutput <- tabOutput[, as.numeric(tableCols)]
  colnames(tabOutput) <- tableNames

  ## return data.frame
  return(tabOutput)

}





## obtain integrated targets to enrich
selectTargets <- function(mirnaObj, miRNA.Direction) {
  
  ## extract integrated targets
  intRes <- integration(mirnaObj)
  if (colnames(intRes)[2] == "Target") {
    targets <- unique(intRes$Target[intRes$microRNA.Direction ==
                                      miRNA.Direction])
  } else if (colnames(intRes)[2] == "mirna.direction") {
    targets <- intRes$DE_targets[intRes$mirna.direction == miRNA.Direction]
    targets <- paste(targets, collapse = "/")
    targets <- strsplit(targets, "/")
    targets <- unlist(targets)
  }
  
  ## return targets
  return(targets)
  
}





## hide output by cat
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}





#' Get the list of supported organisms for a given database
#'
#' This function provides the list of supported organisms for different
#' databases, namely Gene Ontology (GO), Kyoto Encyclopedia of Genes and
#' Genomes (KEGG), MsigDB, WikiPathways, Reactome, Enrichr, Disease Ontology
#' (DO), Network of Cancer Genes (NCG), DisGeNET, and COVID19.
#'
#' @param database The database name. It must be one of: `GO`, `KEGG`, `MsigDB`,
#' `WikiPathways`, `Reactome`, `Enrichr`, `DO`, `NCG`, `DisGeNET`, `COVID19`
#'
#' @returns
#' A `character` vector listing all the supported organisms for the database
#' specified by the user.
#'
#' @examples
#' # get the supported organisms for GO database
#' supportedOrganisms("GO")
#'
#' # get the supported organisms for Reactome
#' supportedOrganisms("Reactome")
#' 
#' @note
#' To perform the functional enrichment of genes, MIRit uses the `geneset` R
#' package to download gene sets from the above mentioned databases.
#'
#' @references
#' Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr R
#' package and web server. BMC Bioinformatics 24, 214 (2023).
#' \url{https://doi.org/10.1186/s12859-023-05342-9}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
supportedOrganisms <- function(database) {
  
  ## check inputs
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "MsigDB", "WikiPathways", "Reactome",
                       "Enrichr", "DO", "NCG", "DisGeNET", "COVID19")) {
    stop(paste("'database' must be one of 'GO', 'KEGG', 'MsigDB',",
               "'WikiPathways', 'Reactome', 'Enrichr', 'DO', 'NCG',",
               "'DisGeNET', 'COVID19'. For additional details,",
               "see ?supportedOrganisms"),
         call. = FALSE)
  }
  
  ## extract supported organisms from species data.frame
  supp <- species[!is.na(species[, database]), "specie"]
  
  ## return supported organisms
  return(supp)
  
}





#' List all the available biological pathways in KEGG, Reactome and
#' WikiPathways
#'
#' This function can be used to retrieve a list of valid biological pathways
#' present in KEGG, Reactome and WikiPathways.
#'
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function
#' @param database The name of the database to use. It must be one of: `KEGG`,
#' `Reactome`, and `WikiPathways`
#'
#' @returns
#' A `character` vector containing the pathway names present in the
#' specified database.
#'
#' @examples
#' # list the mouse pathways present in WikiPathways
#' listPathways("Mus musculus", "WikiPathways")
#'
#' @note
#' This function uses the `graphite` package to retrieve biological pathways
#' from KEGG, Reactome and WikiPathways.
#'
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012),
#' \url{https://doi.org/10.1186/1471-2105-13-20}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
listPathways <- function(organism, database) {
  
  ## check inputs
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("KEGG", "Reactome", "WikiPathways")) {
    stop("Supported databases are: 'KEGG', 'Reactome' and 'WikiPathways'",
         call. = FALSE)
  }
  
  ## check if database is supported for the given specie
  supp <- species[!is.na(species[, paste("graph", database, sep = "_")]),
                  "specie"]
  if (!organism %in% supp) {
    stop(paste("For", database, "database, 'organism' must be one of:",
               paste(supp, collapse = ", ")))
  }
  
  ## set organism name
  organism <- species[species$specie == organism,
                      paste("graph", database, sep = "_")]
  
  ## download pathways from specified database
  pathDb <- graphite::pathways(species = organism,
                               database = tolower(database))
  
  ## return the names of the pathways present in the specified database
  return(names(pathDb))
  
}





#' Create example [`MirnaExperiment`][MirnaExperiment-class] objects
#' 
#' This helper function allows to create a
#' [`MirnaExperiment`][MirnaExperiment-class] containing miRNA and gene
#' expression data deriving from Riesco-Eizaguirre et al (2015).
#'
#' @returns
#' An example `MirnaExperiment` object.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
loadExamples <- function() {
  
  ## return the example object
  return(exampleObject)
  
}

