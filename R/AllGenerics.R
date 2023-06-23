## ==========================================================================
## Generics for MirnaExperiment class
## ==========================================================================



#' Extract differentially expressed miRNAs from a
#' [`MirnaExperiment`][MirnaExperiment-class] object
#'
#' This function is an accessor for the `mirnaDE` slot of
#' [`MirnaExperiment`][MirnaExperiment-class] class, and thus, it can be used
#' to retrieve the results of miRNA differential expression stored in a
#' [`MirnaExperiment`][MirnaExperiment-class] object.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#' @param onlySignificant Logical, if `TRUE` differential expression results
#' will be returned just for statistically significant miRNAs, if `FALSE` the
#' full table of miRNA differential expression will be provided. Default is
#' `TRUE` to only report significant miRNAs
#' @param param Logical, whether to return the complete `list` object with
#' the parameters used, or just the results stored in `data`. Default is FALSE
#' @param returnObject Logical, if `TRUE` this function will return the
#' limma/edgeR/DESeq2 object used for differential expression analysis
#'
#' @returns
#' A `data.frame` with miRNA differential expression, or a `list` object with
#' the parameters used if `param = TRUE`.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # access miRNA differential expression of a MirnaExperiment object
#' sig <- mirnaDE(obj)
#' all <- mirnaDE(obj, onlySignificant = FALSE)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("mirnaDE", function(object,
                               onlySignificant = TRUE,
                               param = FALSE,
                               returnObject = FALSE)
  standardGeneric("mirnaDE"))



#' Extract differentially expressed genes from a
#' [`MirnaExperiment`][MirnaExperiment-class] object
#'
#' This function is an accessor for the `geneDE` slot of
#' [`MirnaExperiment`][MirnaExperiment-class]  class, and thus, it can be used
#' to retrieve the results of gene differential expression stored in a
#' [`MirnaExperiment`][MirnaExperiment-class] object.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#' @param onlySignificant Logical, if `TRUE` differential expression results
#' will be returned just for statistically significant genes, if `FALSE` the
#' full table of gene differential expression will be provided. Default is
#' `TRUE` to only report significant genes
#' @param param Logical, whether to return the complete `list` object with
#' the parameters used, or just the results stored in `data`. Default is FALSE
#' @param returnObject Logical, if `TRUE` this function will return the
#' limma/edgeR/DESeq2 object used for differential expression analysis
#'
#' @returns
#' A `data.frame` with gene differential expression, or a `list` object with
#' the parameters used if `param = TRUE`.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # access gene differential expression of a MirnaExperiment object
#' sig <- geneDE(obj)
#' all <- geneDE(obj, onlySignificant = FALSE)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("geneDE", function(object,
                              onlySignificant = TRUE,
                              param = FALSE,
                              returnObject = FALSE)
  standardGeneric("geneDE"))



#' Get the IDs of statistically differentially expressed miRNAs
#'
#' This function accesses the `significant` features contained in the `mirnaDE`
#' slot of a [`MirnaExperiment`][MirnaExperiment-class] object, and can be used
#' to obtain the IDs of statistically differentially expressed miRNAs.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `character` vector of miRNA IDs (e.g. 'hsa-miR-16-5p',
#' 'hsa-miR-29a-3p'...).
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # extract significant DE-miRNAs
#' sigMirnas <- significantMirnas(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("significantMirnas", function(object)
  standardGeneric("significantMirnas"))



#' Get the IDs of statistically differentially expressed genes
#'
#' This function accesses the `significant` features contained in the `geneDE`
#' slot of a [`MirnaExperiment`][MirnaExperiment-class] object, and can be used
#' to obtain the IDs of statistically differentially expressed genes.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `character` vector of gene symbols (e.g. 'TP53', 'FOXP2', 'TIGAR'
#' 'CASP1'...).
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # extract significant DEGs
#' sigGenes <- significantGenes(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("significantGenes", function(object)
  standardGeneric("significantGenes"))



#' View the relationship between miRNA and gene samples
#'
#' This function allows to access the `pairedSamples` slot of a
#' [`MirnaExperiment`][MirnaExperiment-class] object. The `MirnaExperiment`
#' class is able to contain miRNA and gene expression measurements deriving
#' from the same individuals (paired samples), or from different subjects
#' (unpaired samples).
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `logical` value that is either `TRUE` for paired samples, or `FALSE` for
#' unpaired samples.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # check if an existing MirnaExperiment object derive from paired samples
#' pairedSamples(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("pairedSamples", function(object) standardGeneric("pairedSamples"))



#' Explore miRNA-target pairs
#'
#' This function accesses the `targets` slot of a
#' [`MirnaExperiment`][MirnaExperiment-class]
#' object. After retrieving miRNA targets with the [getTargets()] function,
#' the interactions between miRNAs and target genes are stored in the
#' `targets` slot and can be explored with this function. If `demTarg`
#' parameter is set to TRUE, only targets of differentially expressed miRNAs
#' will be considered; otherwise, all miRNA-target relationships will be shown.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#' @param demTarg Logical, whether to report the targets of differentially
#' expressed miRNAs, or to report the targets of all miRNAs in study.
#' Default is TRUE to just report the target genes of differentially
#' expressed miRNAs
#'
#' @returns
#' A `data.frame` object containing the interactions between miRNAs and target
#' genes, as retrieved with the [getTargets()] function.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # visualize targets
#' targets_df <- mirnaTargets(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("mirnaTargets",
           function(object,
                    demTarg = TRUE) standardGeneric("mirnaTargets"))



#' Explore the results of the integration analysis between miRNAs and genes
#'
#' After performing the integration analysis between miRNA and gene expression
#' values with the [mirnaIntegration()] function, the results are stored
#' in the `integration` slot of a
#' [`MirnaExperiment`][MirnaExperiment-class] object and can be explored with
#' this function.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#' @param param Logical, whether to return the complete `list` object with
#' the parameters used, or just the results stored in `data`. Default is FALSE
#'
#' @returns
#' If `param` is FALSE, then this functions returns a `data.frame` object
#' containing the results of the integration analysis. Otherwise, a `list`
#' object including the parameters used for the analysis will be returned.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform Kendall's correlation analysis with tau > 0.8 and p < 0.05
#' obj <- mirnaIntegration(obj, test = "correlation",
#' corMethod = "kendall", corCutoff = 0.8)
#'
#' # visualize the results of correlation analysis
#' res <- integration(obj)
#' res
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("integration",
           function(object,
                    param = FALSE) standardGeneric("integration"))



setGeneric("mirnaDE<-", function(object, value)
  standardGeneric("mirnaDE<-")
)

setGeneric("geneDE<-", function(object, value)
  standardGeneric("geneDE<-")
)

setGeneric("pairedSamples<-", function(object, value)
  standardGeneric("pairedSamples<-")
)

setGeneric("mirnaTargets<-", function(object, value)
  standardGeneric("mirnaTargets<-")
)

setGeneric("integration<-", function(object, value)
  standardGeneric("integration<-")
)



## ==========================================================================
## Generics for FunctionalEnrichment class
## ==========================================================================



#' Extract results from enrichment objects
#'
#' This function allows to access the `data` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object. This is useful 
#' to take a closer look at all the enriched terms of an enrichment analysis. Nfor objects of classes
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results.
#'
#' @returns
#' A `data.frame` object with the full results of the enrichment analysis.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG database
#' gse <- enrichMirnas(obj, organism = "Homo sapiens",
#' database = "KEGG", method = "GSEA")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' enrichmentDotplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentResults", function(object)
  standardGeneric("enrichmentResults"))

setGeneric("enrichmentDatabase", function(object)
  standardGeneric("enrichmentDatabase"))

setGeneric("enrichmentMethod", function(object)
  standardGeneric("enrichmentMethod"))

setGeneric("geneSet", function(object)
  standardGeneric("geneSet"))

setGeneric("enrichmentMetric", function(object)
  standardGeneric("enrichmentMetric"))

setGeneric("enrichedFeatures", function(object)
  standardGeneric("enrichedFeatures"))

setGeneric("enrichmentResults<-", function(object, value)
  standardGeneric("enrichmentResults<-")
)

setGeneric("enrichmentDatabase<-", function(object, value)
  standardGeneric("enrichmentDatabase<-")
)

