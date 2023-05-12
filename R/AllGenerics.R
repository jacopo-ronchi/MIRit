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
#' @param returnObject Logical, if `TRUE` this function will return the
#' limma/edgeR/DESeq2 object used for differential expression analysis
#'
#' @returns
#' A `data.frame` with miRNA differential expression.
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
#' @param returnObject Logical, if `TRUE` this function will return the
#' limma/edgeR/DESeq2 object used for differential expression analysis
#'
#' @returns
#' A `data.frame` with gene differential expression.
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



#' Explore the targets of differentially expressed microRNAs
#'
#' This function accesses the `mirnaTargets` slot of a
#' [`MirnaExperiment`][MirnaExperiment-class]
#' object. After retrieving miRNA targets with the [getTargets()] function,
#' the interactions between miRNAs and target genes are stored in the
#' `mirnaTargets` slot and can be explored with this function.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `data.frame` object containing the interactions (predicted and/or
#' validated) between miRNAs and target genes, as retrieved with the
#' [getTargets()] function.
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
setGeneric("mirnaTargets", function(object) standardGeneric("mirnaTargets"))



#' Explore the results of the integration analysis between miRNAs and genes
#'
#' After performing the integration analysis between miRNA and gene expression
#' values with the [integrateMirnaTargets()] function, the results are stored
#' in the `mirnaTargetsIntegration` slot of a
#' [`MirnaExperiment`][MirnaExperiment-class] object and can be explored with
#' this function.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `data.frame` object containing the results of the integration analysis.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform Kendall's correlation analysis with tau > 0.8 and p < 0.05
#' obj <- integrateMirnaTargets(obj, test = "correlation",
#' corMethod = "kendall", corCutoff = 0.8)
#'
#' # visualize the results of correlation analysis
#' res <- mirnaTargetsIntegration(obj)
#' res
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("mirnaTargetsIntegration",
           function(object) standardGeneric("mirnaTargetsIntegration"))



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

setGeneric("mirnaTargetsIntegration<-", function(object, value)
  standardGeneric("mirnaTargetsIntegration<-")
)



## ==========================================================================
## Generics for microRNA enrichment classes: MirnaEnrichment & MirnaGsea
## ==========================================================================



#' Create a dotplot for miRNA enrichment analysis
#'
#' This function produces a dotplot to show the results of miRNA enrichment
#' analysis. This function can be used to plot results of both
#' over-representation analysis (ORA), carried out with [enrichMirnas()], and
#' gene set enrichment analysis (GSEA), performed with [gseaMirnas()].
#'
#' @param object An object of class [`MirnaEnrichment`][MirnaEnrichment-class]
#' or [`MirnaGsea`][MirnaGsea-class] containing miRNA enrichment results
#' @param showTerms It is the number of top enriched terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param splitDir Logical, if `TRUE` the resulting plot will be divided in
#' two columns on the basis of enrichment direction (enrichment and depletion).
#' Default is `FALSE`
#' @param ordBy The parameter used to set the x-axis scale. It must be one of
#' `fold` (default), `P.adjusted`, `P.value` and `Observed`
#' @param sizeBy The parameter used to set the size scale. It must be one of
#' `fold`, `P.adjusted`, `P.value` and `Observed` (default)
#' @param colBy The parameter used to set the color scale. It must be one of
#' `fold`, `P.adjusted` (default), `P.value` and `Observed`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @note
#' For objects of class [`MirnaEnrichment`][MirnaEnrichment-class], the `fold`
#' value in `ordBy`, `sizeBy` and `colBy` refers to fold enrichment, which
#' corresponds to the ratio between observed and expected hits for a given
#' category. This means that the higher the fold enrichment is, the more
#' unlikely it is that the enrichment was reported by chance. Instead, for
#' objects of class [`MirnaGsea`][MirnaGsea-class] the `fold` value refers to
#' the average absolute fold change of observed miRNAs in a given category.
#'
#' @returns
#' A `ggplot` graph with a dotplot of enrichment results.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform a GSEA with the miEAA category GO Annotations
#' gse <- gseaMirnas(obj, organism = "Homo sapiens",
#' category = "GO_Annotations_mature")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' mirnaDotplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("mirnaDotplot",
           function(object,
                    showTerms = 10,
                    splitDir = FALSE,
                    ordBy = "fold",
                    sizeBy = "Observed",
                    colBy = "P.adjusted",
                    title = NULL) standardGeneric("mirnaDotplot"))



#' Extract results from enrichment objects
#'
#' This function allows to access the `result` slot of
#' [`MirnaEnrichment`][MirnaEnrichment-class] and [`MirnaGsea`][MirnaGsea-class]
#' objects. This is useful to take a closer look at all the
#' enriched terms of an enrichment analysis. Note that this function can also
#' extract the enrichment results for objects of classes
#' [`enrichResult-class`][DOSE::enrichResult-class] and
#' [`gseaResult-class`][DOSE::gseaResult-class].
#'
#' @param object An object of class [`MirnaEnrichment`][MirnaEnrichment-class]
#' or [`MirnaGsea`][MirnaGsea-class] containing miRNA enrichment results.
#' Alternatively, object of class
#' [`enrichResult-class`][DOSE::enrichResult-class] or
#' [`gseaResult-class`][DOSE::gseaResult-class] are also accepted
#'
#' @returns
#' A `data.frame` object with the full results of the enrichment analysis.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform a GSEA with the miEAA category GO Annotations
#' gse <- gseaMirnas(obj, organism = "Homo sapiens",
#' category = "GO_Annotations_mature")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' mirnaDotplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentResults", function(object)
  standardGeneric("enrichmentResults"))


## Make 'enrichmentResults' work also for enrichResult and gseaResult
#' @rdname enrichmentResults
#' @export
setMethod("enrichmentResults", "enrichResult", function(object) {
  object@result
})

#' @rdname enrichmentResults
#' @export
setMethod("enrichmentResults", "gseaResult", function(object) {
  object@result
})



setGeneric("enrichmentDatabase", function(object)
  standardGeneric("enrichmentDatabase"))

setGeneric("lfcEnrichment", function(object)
  standardGeneric("lfcEnrichment"))

setGeneric("mirnaIdEnrichment", function(object)
  standardGeneric("mirnaIdEnrichment"))

setGeneric("enrichmentResults<-", function(object, value)
  standardGeneric("enrichmentResults<-")
)

setGeneric("enrichmentDatabase<-", function(object, value)
  standardGeneric("enrichmentDatabase<-")
)


