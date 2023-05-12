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
      nrow(mirnaTargets(object)), " rows and ", ncol(mirnaTargets(object)),
      " columns\n",
      "\t- miRNA - gene integrative analysis: ",
      class(mirnaTargetsIntegration(object)), " with ",
      nrow(mirnaTargetsIntegration(object)), " rows and ",
      ncol(mirnaTargetsIntegration(object)), " columns\n",
      "\nMicroRNA and gene expression data derive from: ",
      ifelse(pairedSamples(object) == TRUE, "paired", "unpaired"),
      " samples\n\n", sep = "")
})



## ==========================================================================
## Show method for MirnaEnrichment class
## ==========================================================================


#' @describeIn MirnaEnrichment-class Show method for objects of
#' class MirnaEnrichment
#' @export
setMethod("show", "MirnaEnrichment", function(object) {
  cat("Object of class MirnaEnrichment containing:\n\n",
      "\t- over-representation analysis results: ",
      class(enrichmentResults(object)), " with ",
      nrow(enrichmentResults(object)), " rows and ",
      ncol(enrichmentResults(object)), " columns\n",
      "\t- miRNAs used for the enrichment: ", class(object@gene), " of length ",
      length(object@gene), "\n",
      "\t- p-value cutoff used: ", object@pvalueCutoff, "\n",
      "\t- p-value adjustment method: ", object@pAdjustMethod, "\n",
      "\t- miEAA category: ", enrichmentDatabase(object), "\n",
      "\t- organism: ", object@organism, "\n",
      "\nResults can be accessed with 'enrichmentResults()' while ",
      "'mirnaDotplot()' can be used for visualization.\n\n", sep = "")
})


## ==========================================================================
## Show method for MirnaGsea class
## ==========================================================================


#' @describeIn MirnaGsea-class Show method for objects of class MirnaGsea
#' @export
setMethod("show", "MirnaGsea", function(object) {
  cat("Object of class MirnaGsea containing:\n\n",
      "\t- GSEA analysis results: ",
      class(enrichmentResults(object)), " with ",
      nrow(enrichmentResults(object)), " rows and ",
      ncol(enrichmentResults(object)), " columns\n",
      "\t- miRNAs used for the enrichment: ", class(object@gene), " of length ",
      length(object@gene), "\n",
      "\t- miRNA log foldchanges: ", class(lfcEnrichment(object)),
      " of length ", length(lfcEnrichment(object)), "\n",
      "\t- p-value cutoff used: ", object@pvalueCutoff, "\n",
      "\t- p-value adjustment method: ", object@pAdjustMethod, "\n",
      "\t- miEAA category: ", enrichmentDatabase(object), "\n",
      "\t- organism: ", object@organism, "\n",
      "\nResults can be accessed with 'enrichmentResults()' while ",
      "'mirnaDotplot()' and 'mirnaRidgeplot()' can be used for ",
      "visualization.\n\n", sep = "")
})

