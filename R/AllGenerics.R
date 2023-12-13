## ==========================================================================
## Generics for MirnaExperiment class
## ==========================================================================

#' Extract differentially expressed miRNAs and genes from a
#' [`MirnaExperiment`][MirnaExperiment-class] object
#'
#' The `mirnaDE()` and `geneDE()` are two accessor functions for the `mirnaDE`
#' and `geneDE` slots of [`MirnaExperiment`][MirnaExperiment-class] class,
#' respectively. Thus, they can be used to explore the results of miRNA and
#' gene differential expression analysis stored in a
#' [`MirnaExperiment`][MirnaExperiment-class] object.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#' @param onlySignificant Logical, if `TRUE` differential expression results
#' will be returned just for statistically significant miRNAs/genes, if `FALSE`
#' the full table of miRNA/gene differential expression will be provided.
#' Default is `TRUE` to only report significant miRNAs/genes
#' @param param Logical, whether to return the complete `list` object with
#' the parameters used, or just the results stored in `data`. Default is FALSE
#' @param returnObject Logical, if `TRUE` this function will return the
#' limma/edgeR/DESeq2 object used for differential expression analysis
#'
#' @returns
#' A `data.frame` with miRNA/gene differential expression, or a `list` object
#' with the parameters used if `param = TRUE`.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # access miRNA differential expression of a MirnaExperiment object
#' sig <- mirnaDE(obj)
#' all <- mirnaDE(obj, onlySignificant = FALSE)
#'
#' # access gene differential expression of a MirnaExperiment object
#' sig <- geneDE(obj)
#' all <- geneDE(obj, onlySignificant = FALSE)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name deAccessors
NULL



#' @describeIn deAccessors Extract miRNA differential expression results
#' @export
setGeneric("mirnaDE", function(object,
    onlySignificant = TRUE,
    param = FALSE,
    returnObject = FALSE) {
    standardGeneric("mirnaDE")
})



#' @describeIn deAccessors Extract gene differential expression results
#' @export
setGeneric("geneDE", function(object,
    onlySignificant = TRUE,
    param = FALSE,
    returnObject = FALSE) {
    standardGeneric("geneDE")
})



#' Get the IDs of statistically differentially expressed miRNAs/genes
#'
#' The `significantMirnas()` and `significantGenes()` functions access the
#' `significant` features contained in the `mirnaDE` or `geneDE` slots
#' of a [`MirnaExperiment`][MirnaExperiment-class] object, and can be used
#' to obtain the IDs of statistically differentially expressed miRNAs and genes.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
#'
#' @returns
#' A `character` vector of miRNA IDs (e.g. 'hsa-miR-16-5p',
#' hsa-miR-29a-3p'...), or a`character` vector of gene symbols (e.g. 'TP53',
#' 'FOXP2', 'TIGAR', CASP1'...).
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # extract significant DE-miRNAs
#' sigMirnas <- significantMirnas(obj)
#'
#' # extract significant DEGs
#' sigGenes <- significantGenes(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name significantAccessors
NULL



#' @describeIn significantAccessors Get the IDs of differentially expressed
#' miRNAs
#' @export
setGeneric("significantMirnas", function(object) {
    standardGeneric("significantMirnas")
})



#' @describeIn significantAccessors Get the IDs of differentially expressed
#' genes
#' @export
setGeneric("significantGenes", function(object) {
    standardGeneric("significantGenes")
})



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
#' `targets` slot and can be explored with this function.
#'
#' @param object A [`MirnaExperiment`][MirnaExperiment-class] object containing
#' miRNA and gene data
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
setGeneric(
    "mirnaTargets",
    function(object) standardGeneric("mirnaTargets")
)



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
#' obj <- mirnaIntegration(obj,
#'     test = "correlation",
#'     corMethod = "kendall", corCutoff = 0.8
#' )
#'
#' # visualize the results of correlation analysis
#' res <- integration(obj)
#' res
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric(
    "integration",
    function(object,
    param = FALSE) {
        standardGeneric("integration")
    }
)



setGeneric("mirnaDE<-", function(object, value) {
    standardGeneric("mirnaDE<-")
})

setGeneric("geneDE<-", function(object, value) {
    standardGeneric("geneDE<-")
})

setGeneric("pairedSamples<-", function(object, value) {
    standardGeneric("pairedSamples<-")
})

setGeneric("mirnaTargets<-", function(object, value) {
    standardGeneric("mirnaTargets<-")
})

setGeneric("integration<-", function(object, value) {
    standardGeneric("integration<-")
})



## ==========================================================================
## Generics for FunctionalEnrichment class
## ==========================================================================

#' Access the results of functional enrichment analyses
#'
#' This function accesses the `data` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `data.frame` with enrichment results.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `data.frame` object containing the results of functional enrichment
#' analyses, as returned by the [enrichGenes()] function.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract results
#' de_df <- enrichmentResults(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentResults", function(object) {
    standardGeneric("enrichmentResults")
})



#' Access the database used for functional enrichment analyses
#'
#' This function accesses the `database` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a the name of the database used by the [enrichGenes()] function to perform
#' the enrichment analysis.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `character` containing the name of the database, such as `KEGG`.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # see the database
#' enrichmentDatabase(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentDatabase", function(object) {
    standardGeneric("enrichmentDatabase")
})



#' Access the method used for functional enrichment analyses
#'
#' This function accesses the `method` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a the name of the enrichment strategy used by the [enrichGenes()] function
#' to perform the enrichment analysis.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `character` containing the enrichment method, such as `GSEA`.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # see the method
#' enrichmentMethod(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentMethod", function(object) {
    standardGeneric("enrichmentMethod")
})



#' Extract the gene-sets used for functional enrichment analyses
#'
#' This function accesses the `geneSet` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `list` with the collection of genes used for the enrichment.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `list` containing the gene-sets.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract the gene-sets
#' gs <- geneSet(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("geneSet", function(object) {
    standardGeneric("geneSet")
})



#' Extract the GSEA ranking metric used for functional enrichment analyses
#'
#' This function accesses the `statistic` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `numeric` vector with the metric used to rank genes in GSEA.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `numeric` vector containing the ranking metric.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract the ranking metric
#' rmet <- enrichmentMetric(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichmentMetric", function(object) {
    standardGeneric("enrichmentMetric")
})



#' Extract the names of the pre-ranked features in a GSEA experiment
#'
#' This function accesses the `features` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `character` vector with the names of the features considered in GSEA in
#' the same order as the ranking metric.
#'
#' @param object An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#'
#' @returns
#' A `character` vector with the names of the genes ordered based on the
#' ranking metric.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract the ranking metric
#' rmet <- enrichmentMetric(obj)
#'
#' ## extract the corresponding names
#' rnames <- enrichedFeatures(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("enrichedFeatures", function(object) {
    standardGeneric("enrichedFeatures")
})



setGeneric("enrichmentResults<-", function(object, value) {
    standardGeneric("enrichmentResults<-")
})

setGeneric("enrichmentDatabase<-", function(object, value) {
    standardGeneric("enrichmentDatabase<-")
})



## ==========================================================================
## Generics for IntegrativePathwayAnalysis class
## ==========================================================================

#' Access the results of integrative miRNA-mRNA pathway analyses
#'
#' This function accesses the `data` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `data.frame` with the results of an integrative topological analysis
#' carried out through the [topologicalAnalysis()] function.
#'
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] containing
#' the results of a miRNA-mRNA pathway analysis
#'
#' @returns
#' A `data.frame` object containing the results of the topological analysis,
#' as returned by the [topologicalAnalysis()] function.
#'
#' @examples
#' # load the example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#'
#' # extract results
#' taipaRes <- integratedPathways(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("integratedPathways", function(object) {
    standardGeneric("integratedPathways")
})



#' Extract the database used for integrative miRNA-mRNA pathway analyses
#'
#' This function accesses the `database` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' the name of the database used by the [topologicalAnalysis()] function to
#' perform the integrative topological analysis.
#'
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] containing
#' the results of a miRNA-mRNA pathway analysis
#'
#' @returns
#' A `character` object with the name of the database used by the
#' [topologicalAnalysis()] function, such as `KEGG`.
#'
#' @examples
#' # load the example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#'
#' # see the database
#' integrationDatabase(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("integrationDatabase", function(object) {
    standardGeneric("integrationDatabase")
})



#' Access the miRNA-augmented pathways that were used during TAIPA
#'
#' This function accesses the `pathways` slot of a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] object and returns
#' a `list` object with the augmented pathways that were considered by the
#' [topologicalAnalysis()] function to perform the integrative analysis.
#'
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] containing
#' the results of a miRNA-mRNA pathway analysis
#'
#' @returns
#' A `list` object with the miRNA-augmented biological pathways.
#'
#' @examples
#' # load the example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#'
#' # extract the pathways
#' ps <- augmentedPathways(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
setGeneric("augmentedPathways", function(object) {
    standardGeneric("augmentedPathways")
})



setGeneric("integratedPathways<-", function(object, value) {
    standardGeneric("integratedPathways<-")
})
