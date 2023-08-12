## ==========================================================================
## Show method for MirnaExperiment class
## ==========================================================================


#' @describeIn MirnaExperiment-class Show method for objects of class
#' MirnaExperiment
#' @export
#' @importMethodsFrom MultiAssayExperiment show
setMethod("show", "MirnaExperiment", function(object) {
  cat("An object of class MirnaExperiment, which extends MultiAssayExperiment",
      "class and contains:\n\n",
      "\t- microRNA expression values: ",
      class(MultiAssayExperiment::experiments(object)[["microRNA"]]), " with ",
      nrow(MultiAssayExperiment::experiments(object)[["microRNA"]]),
      " rows and ",
      ncol(MultiAssayExperiment::experiments(object)[["microRNA"]]),
      " columns\n", "\t- gene expression values: ",
      class(MultiAssayExperiment::experiments(object)[["genes"]]), " with ",
      nrow(MultiAssayExperiment::experiments(object)[["genes"]]), " rows and ",
      ncol(MultiAssayExperiment::experiments(object)[["genes"]]), " columns\n",
      "\t- samples metadata: ", class(MultiAssayExperiment::colData(object)),
      " with ", nrow(MultiAssayExperiment::colData(object)), " rows and ",
      ncol(MultiAssayExperiment::colData(object))," columns\n",
      "\t- microRNA differential expression: ",
      class(mirnaDE(object, onlySignificant = FALSE)),
      " with ", nrow(mirnaDE(object, onlySignificant = FALSE)),
      " rows and ", ncol(mirnaDE(object, onlySignificant = FALSE)),
      " columns\n",
      "\t- significant DE-miRNAs: ", class(significantMirnas(object)), " with ",
      length(significantMirnas(object)), " miRNA IDs\n",
      "\t- gene differential expression: ",
      class(geneDE(object, onlySignificant = FALSE)), " with ",
      nrow(geneDE(object, onlySignificant = FALSE)), " rows and ",
      ncol(geneDE(object, onlySignificant = FALSE)), " columns\n",
      "\t- significant genes: ", class(significantGenes(object)), " with ",
      length(significantGenes(object)), " gene IDs\n",
      "\t- microRNA targets: ", class(mirnaTargets(object)), " with ",
      nrow(mirnaTargets(object)), " rows and ",
      ncol(mirnaTargets(object)), " columns\n",
      "\t- miRNA - gene integrative analysis: ",
      class(integration(object)), " with ",
      nrow(integration(object)), " rows and ",
      ncol(integration(object)), " columns\n",
      "\nMicroRNA and gene expression data derive from: ",
      ifelse(pairedSamples(object) == TRUE, "paired", "unpaired"),
      " samples\n\n", sep = "")
})



## ==========================================================================
## Show method for FunctionalEnrichment class
## ==========================================================================


#' @describeIn FunctionalEnrichment-class Show method for objects of
#' class FunctionalEnrichment
#' @export
setMethod("show", "FunctionalEnrichment", function(object) {
  cat("Object of class FunctionalEnrichment containing:\n\n",
      "\t- over-representation analysis results: ",
      class(enrichmentResults(object)), " with ",
      nrow(enrichmentResults(object)), " rows and ",
      ncol(enrichmentResults(object)), " columns\n",
      "\t- functional enrichment analysis: ", object@method, "\n",
      "\t- organism: ", object@organism, "\n",
      "\t- gene sets database: ", enrichmentDatabase(object), "\n",
      "\t- p-value cutoff used: ", object@pCutoff, "\n",
      "\t- p-value adjustment method: ", object@pAdjustment, "\n",
      "\t- features used for the enrichment: ", class(object@features),
      " of length ", length(object@features), "\n",
      "\t- statistic used (only for GSEA): ", class(object@statistic),
      " of length ", length(object@statistic), "\n",
      "\t- background universe for the enrichment: ", class(object@universe),
      " of length ", length(object@universe), "\n",
      "\nResults can be accessed with 'enrichmentResults()' while ",
      "'enrichmentDotplot()' can be used for visualization.\n\n", sep = "")
})

