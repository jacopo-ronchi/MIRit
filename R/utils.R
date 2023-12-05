## helper function to identify column names
identifyColNames <- function(tabOutput, tabID = "") {
    # define accepted column names
    idAccepted <- "ID|Symbol|Gene_Symbol|Mirna|mir|Gene|gene.symbol|Gene.symbol"
    fcAccepted <- "logFC|log2FoldChange|FC|lFC"
    exprAccepted <- "AveExpr|baseMean|logCPM"
    pvalAccepted <- "P.Value|pvalue|PValue|Pvalue"
    fdrAccepted <- "adj.P.Val|padj|FDR|fdr|adj|adj.p|adjp"
    acceptedNames <- c(
        idAccepted, fcAccepted, exprAccepted,
        pvalAccepted, fdrAccepted
    )
    
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
        stop("MIRit is unable to automatically find columns ",
             "relative to ",
             paste(tableNames[which(lengths(tableCols) == 0)],
                   collapse = ", "
             ),
             " in the ", tabID, " differential expression data.frame! ",
             "Please change the name of these columns to one of: ",
             paste(acceptedNames[which(lengths(tableCols) == 0)],
                   collapse = ", "
             ),
             call. = FALSE
        )
    } else if (any(lengths(tableCols) > 1)) {
        stop("More than one column can be interpreted as ",
             paste(tableNames[which(lengths(tableCols) > 1)],
                   collapse = ", "
             ),
             " in the ", tabID, " differential expression data.frame! ",
             "Please rename these columns to have unambiguous column names! ",
             "See ?MirnaExperiment for further details",
             call. = FALSE
        )
    }
    
    ## create a data.frame with desired columns
    tabOutput <- tabOutput[, as.numeric(tableCols)]
    colnames(tabOutput) <- tableNames
    
    ## return data.frame
    return(tabOutput)
}





## helper function to obtain integrated targets to enrich
selectTargets <- function(mirnaObj, miRNA.Direction = NULL, pairs = FALSE) {
    ## extract integration results
    intRes <- integration(mirnaObj)
    
    ## determine the integration method
    method <- integration(mirnaObj, param = TRUE)$method
    
    ## extract integrated miRNA-target pairs
    if (grepl("correlation", method) == TRUE) {
        intPairs <- intRes[, c(1, 2, 3)]
    } else {
        intRes <- intRes[intRes$DE.targets != "", ]
        intPairs <- do.call(
            rbind,
            apply(intRes, 1, function(x) {
                data.frame(
                    microRNA = x["microRNA"],
                    Target = strsplit(x["DE.targets"], "/"),
                    microRNA.Direction = x["mirna.direction"],
                    row.names = NULL
                )
            })
        )
        colnames(intPairs)[2] <- "Target"
    }
    
    ## select miRNA-target pairs that vary in a specific direction
    if (!is.null(miRNA.Direction)) {
        if (miRNA.Direction == "downregulated") {
            intPairs <- intPairs[
                intPairs$microRNA.Direction == "Down" |
                    intPairs$microRNA.Direction == "downregulated", ]
        } else {
            intPairs <- intPairs[
                intPairs$microRNA.Direction == "Up" |
                    intPairs$microRNA.Direction == "upregulated", ]
        }
    }
    
    ## retain just target names
    if (pairs == FALSE) {
        intPairs <- unique(intPairs$Target)
    }
    
    ## return integrated pairs
    return(intPairs)
}





## helper function for hiding output by cat
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}





## helper function for loading cache
.get_cache <- function() {
    cache <- tools::R_user_dir("MIRit", which = "cache")
    BiocFileCache::BiocFileCache(cache, ask = FALSE)
}





## helper function for checking the validity of color names
areColors <- function(x) {
    vapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)),
                 error = function(e) FALSE
        )
    }, FUN.VALUE = logical(1))
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
        !database %in% c(
            "GO", "KEGG", "MsigDB", "WikiPathways", "Reactome",
            "Enrichr", "DO", "NCG", "DisGeNET", "COVID19"
        )) {
        stop("'database' must be one of 'GO', 'KEGG', 'MsigDB', ",
             "'WikiPathways', 'Reactome', 'Enrichr', 'DO', 'NCG', ",
             "'DisGeNET', 'COVID19'. For additional details, ",
             "see ?supportedOrganisms",
             call. = FALSE
        )
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
             call. = FALSE
        )
    }
    
    ## check if database is supported for the given specie
    supp <- species[
        !is.na(species[, paste("graph", database, sep = "_")]),
        "specie"
    ]
    if (!organism %in% supp) {
        stop("For ", database, " database, 'organism' must be one of: ",
             paste(supp, collapse = ", "),
             call. = FALSE
        )
    }
    
    ## set organism name
    organism <- species[
        species$specie == organism,
        paste("graph", database, sep = "_")
    ]
    
    ## download pathways from specified database
    pathDb <- graphite::pathways(
        species = organism,
        database = tolower(database)
    )
    
    ## return the names of the pathways present in the specified database
    return(names(pathDb))
}





## helper function for pathway node conversion
## this function has been adapted from the graphite package to run via
## BiocParallel parallelization
convertNodes <- function(x, db) {
    ## define conversion mapping
    mapping <- list(to = "SYMBOL", db = db)
    
    ## convert node identifiers
    x@protEdges <- graphite:::convertEdges(x@protEdges, mapping)
    x@protPropEdges <- graphite:::convertEdges(x@protPropEdges, mapping)
    x@metabolEdges <- graphite:::convertEdges(x@metabolEdges, mapping)
    x@metabolPropEdges <- graphite:::convertEdges(x@metabolPropEdges, mapping)
    x@mixedEdges <- graphite:::convertEdges(x@mixedEdges, mapping)
    
    if (nrow(x@protEdges) + nrow(x@protPropEdges) + nrow(x@metabolEdges) +
        nrow(x@metabolPropEdges) + nrow(x@mixedEdges) == 0) {
        warning("the conversion lost all edges of pathway \"", x@title, "\"")
    }
    
    ## return the pathway with converted nodes
    return(x)
}





#' Load example MIRit objects
#'
#' This helper function allows to create a
#' [`MirnaExperiment`][MirnaExperiment-class] object containing miRNA and gene
#' expression data deriving from Riesco-Eizaguirre et al (2015), an
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] object
#' containing TAIPA results for the same dataset, or a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] with example GSEA
#' enrichment results.
#'
#' @param class It must be `MirnaExperiment` (default) to load an example
#' object of class [`MirnaExperiment`][MirnaExperiment-class],
#' `IntegrativePathwayAnalysis`, to load an example object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class], or
#' `FunctionalEnrichment`, to load an example object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class].
#'
#' @returns
#' An example `MirnaExperiment` object, an `IntegrativePathwayAnalysis`
#' object, or a `FunctionalEnrichment` object.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # load example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#'
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
loadExamples <- function(class = "MirnaExperiment") {
    ## check input
    if (!class %in% c(
        "MirnaExperiment", "IntegrativePathwayAnalysis",
        "FunctionalEnrichment"
    )) {
        stop("'class' must be one of 'MirnaExperiment', ",
             "'IntegrativePathwayAnalysis', and 'FunctionalEnrichment'",
             call. = FALSE
        )
    }
    
    ## return the example object
    if (class == "MirnaExperiment") {
        return(exampleObject)
    } else if (class == "IntegrativePathwayAnalysis") {
        return(exampleTAIPA)
    } else if (class == "FunctionalEnrichment") {
        return(exampleGSEA)
    }
}
