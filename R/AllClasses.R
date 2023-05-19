## Load generics before classes
#' @include AllGenerics.R
NULL

## ==========================================================================
## MirnaExperiment Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The 'MirnaExperiment' class
#'
#' @description
#' This class extends the
#' [`MultiAssayExperiment`][MultiAssayExperiment::MultiAssayExperiment-class]
#' from the homonym package to provide the flexibility of handling genomic data
#' of both microRNAs and their targets, allowing to store information about
#' microRNA and gene expression, differential expression results, microRNA
#' targets and miRNA-gene integration analysis.
#'
#' This class can be used to manage genomic data deriving from different
#' sources, like RNA-Seq, microarrays and mass spectrometry. Moreover, microRNA
#' and gene expression levels may originate from the same individuals
#' (paired samples) or from different subjects (unpaired samples).
#'
#' @slot ExperimentList An `ExperimentList` class object for each assay dataset
#' @slot colData A `DataFrame` of all clinical/specimen data available across
#' experiments
#' @slot sampleMap A `DataFrame` of translatable identifiers of samples and
#' participants
#' @slot metadata Additional data describing the object
#' @slot drops A metadata `list` of dropped information
#' @slot mirnaDE A `list` object containing the results of miRNA differential
#' expression
#' @slot geneDE A `list` object containing the results of gene differential
#' expression
#' @slot pairedSamples A `logical` parameter that specifies whether miRNA and
#' gene expression measurements derive from the same individuals (`TRUE`) or
#' from different subjects (`FALSE`)
#' @slot targets A `data.frame` object containing miRNA-target pairs. This
#' slot is commonly populated by the [getTargets()] function
#' @slot mirnaTargetsIntegration A `data.frame` object containing the results
#' of the integration analysis between miRNA and gene expression values. This
#' slot is commonly populated by the [integrateMirnaTargets()] function
#' 
#' @section ExperimentList:
#' 
#' The [`ExperimentList`][MultiAssayExperiment::ExperimentList] slot is designed
#' to contain results from each experiment/assay. In this case, it holds miRNA
#' and gene expression matrices. It contains a
#' [`SimpleList-class`][S4Vectors::SimpleList-class].
#' 
#' @section colData:
#' 
#' The `colData` slot is a collection of primary specimen data valid across all
#' experiments. This slot is strictly of class [`DataFrame`] but arguments for
#' the constructor function allow arguments to be of class `data.frame` and
#' subsequently coerced.
#' 
#' @section sampleMap:
#' 
#' The `sampleMap` contains a [`DataFrame`] of translatable identifiers of
#' samples and participants or biological units. The standard column names of
#' the sampleMap are "assay", "primary", and "colname". Note that the "assay"
#' column is a factor corresponding to the names of each experiment in the
#' [`ExperimentList`][MultiAssayExperiment::ExperimentList]. In the case where
#' these names do not match between the sampleMap and the experiments, the
#' documented experiments in the `sampleMap` take precedence and experiments
#' are dropped by the harmonization procedure. The constructor function will
#' generate a `sampleMap` in the case where it is not provided and this method
#' may be a 'safer' alternative for creating the `MultiAssayExperiment` (so
#' long as the rownames are identical in the `colData`, if provided). An empty
#' `sampleMap` may produce empty experiments if the levels of the "assay"
#' factor in the `sampleMap` do not match the names in the
#' [`ExperimentList`][MultiAssayExperiment::ExperimentList].
#'
#' @section mirnaDE and geneDE:
#' 
#' `mirnaDE` and `geneDE` consist of two `list` objects storing information
#' regarding miRNA and gene differential expression, including:
#' 
#' * `data`, which contains differential expression results in a `data.frame`
#' with five columns:
#' 
#'    + `ID`: indicates the name of the miRNA/gene;
#'    + `logFC`: indicates the fold change of each feature in logarithmic scale;
#'    + `AveExpr`: represents the average expression of each miRNA/gene;
#'    + `P.Value`: indicates the resulting p-value;
#'    + `adj.P.Val`: contains the p-values adjusted for multiple testing.
#' 
#' * `significant`, which is a `character` vector containing the names of
#' significantly differentially expressed miRNAs/genes that passed the
#' thresholds;
#' * `method`, which specifies the procedure used to determine differentially
#' expressed miRNAs/gens (eg. "limma-voom", "edgeR", "DESeq2", "limma");
#' * `pCutoff`, which indicates the p-value cutoff used for DE analysis;
#' * `pAdjustment`, the approach used for multiple testing correction;
#' * `logFC`, which states the log2 Fold Change cutoff used for DE analysis;
#' * `deObject`, an object deriving from limma/edgeR/DESeq2, that holds
#' additional information regarding data processing.
#' 
#' MiRNA differential expression results can be accessed through the [mirnaDE()]
#' function, for additional details see `?mirnaDE`. Similarly, gene
#' differential expression results can be accessed through the [geneDE()]
#' function, for additional details see `?geneDE`.
#'
#' @section pairedSamples:
#'
#' As already mentioned, `pairedSamples` must be `TRUE` when miRNA and gene
#' expression derive from the same subjects, while it must `FALSE` if this is
#' not the case.
#'
#' @section targets:
#'
#' `targets` is a `data.frame` with just two columns:
#' * `mature_mirna_id`, which contains miRNA names; and
#' * `target_symbol`, which indicates the target gene for the corresponding
#' miRNA.
#'
#' @section mirnaTargetsIntegration:
#'
#' Lastly, `mirnaTargetsIntegration` slot contains a `data.frame` that differs
#' on the basis of the integration strategy used. For the one-sided Fisher's
#' exact test integration, this `data.frame` has seven columns:
#' * `microRNA`: the miRNA ID;
#' * `direction`: the fold change direction of the DE-miRNA (`upregulated` or
#' `downregulated`);
#' * `n_DE_targets`: represents the number of differentially expressed targets;
#' * `n_NON_DE_targets`: represents the number of non differentially expressed
#' targets;
#' * `Fisher.P.Val`: indicates the p-value resulting from the one-sided
#' Fisher's exact test;
#' * `Fisher.Adjusted.P.Val`: contains the Fisher's exact test p-values
#' corrected for multiple testing;
#' * `DE_targets`: contains the list of differentially expressed targets whose
#' expression is negatively associated with miRNA expression.
#' Instead, when a correlation analysis is performed, `mirnaTargetsIntegration`
#' has seven columns:
#' * `microRNA`: the miRNA ID;
#' * `Target`: the correlated target gene;
#' * `microRNA.Direction`: the fold change direction of the DE-miRNA;
#' * `Pearson/Spearman/Kendall.Coeff`: the value of the correlation coefficient
#' used;
#' * `Corr.P.Value`: the p-value resulting from the correlation analysis;
#' * `Corr.Adjusted.P.Val`: contains the correlation p-values corrected for
#' multiple testing.
#'
#' @note
#' To create a [`MirnaExperiment`][MirnaExperiment-class] object, you can use
#' the [MirnaExperiment()] constructor function, which allows to easily build
#' and verify a valid object starting from miRNA and gene expression matrices.
#' 
#' @param object An object of class [`MirnaExperiment`][MirnaExperiment-class]
#'
#' @name MirnaExperiment-class
#'
#' @references
#' Marcel Ramos et al. Software For The Integration Of Multiomics Experiments
#' In Bioconductor. Cancer Research, 2017 November 1; 77(21); e39-42. DOI:
#' \url{10.1158/0008-5472.CAN-17-0344}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @seealso
#' See [`MultiAssayExperiment-class`][MultiAssayExperiment::MultiAssayExperiment-class]
#' for additional information.
#'
#' @docType class
#' @export
#' @import methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setClass("MirnaExperiment",
         contains = "MultiAssayExperiment",
         slots = representation(
           mirnaDE = "list",
           geneDE = "list",
           pairedSamples = "logical",
           targets = "data.frame",
           mirnaTargetsIntegration = "data.frame"))


## -----------------
## Initialize method
## -----------------


setMethod("initialize",
          signature(.Object = "MirnaExperiment"),
          function(.Object, ...) {
            requiredObjects <- c("data", "significant", "method", "pCutoff",
                                 "pAdjustment", "logFC", "deObject")
            deList <- list(data = data.frame(),
                           significant = character(),
                           method = character(),
                           pCutoff = numeric(),
                           pAdjustment = character(),
                           logFC = numeric(),
                           deObject = NULL)
            
            .Object <- callNextMethod(.Object,
                                      ...,
                                      mirnaDE = deList,
                                      geneDE = deList)
            .Object
          })


## --------
## Validity
## --------

setValidity("MirnaExperiment", function(object) {
  
  if (!is.list(object@mirnaDE)) {
    return(paste("'mirnaDE' slot must be a list object with miRNA",
                 "differential expression results. Please see",
                 "?MirnaExperiment-class"))
  } else if (!is.list(object@geneDE)) {
    return(paste("'geneDE' slot must be a list object with gene",
                 "differential expression results. Please see",
                 "?MirnaExperiment-class"))
  } else if (!identical(sort(names(object@mirnaDE)),
                        sort(c("data", "significant", "method", "pCutoff",
                               "pAdjustment", "logFC", "deObject"))) |
             !identical(sort(names(object@geneDE)),
                        sort(c("data", "significant", "method", "pCutoff",
                               "pAdjustment", "logFC", "deObject")))) {
    return(paste("'mirnaDE' and 'geneDE' slots must be list objects",
                 "containing: 'data', 'significant', 'method', 'pCutoff'",
                 "'pAdjustment', 'logFC', and 'deObject'.",
                 "Please see ?MirnaExperiment-class"))
  } else if (!is.data.frame(mirnaDE(object))) {
    return(paste("'data' within mirnaDE slot must be a data.frame object with",
                 "miRNA differential expression results. Please see",
                 "?MirnaExperiment-class"))
  } else if (!is.data.frame(geneDE(object))) {
    return(paste("'data' within geneDE slot must be a data.frame object with",
                 "gene differential expression results. Please see",
                 "?MirnaExperiment-class"))
  } else if (!is.character(significantMirnas(object))  |
             !all(significantMirnas(object) %in%
                  mirnaDE(object, onlySignificant = FALSE)$ID)) {
    return(paste("'significant' object within 'mirnaDE' slot must be a",
                 "character with IDs of statiscally significantly",
                 "differentially expressed miRNAs."))
  } else if (!is.character(significantGenes(object))  |
             !all(significantGenes(object)
                  %in% geneDE(object, onlySignificant = FALSE)$ID)) {
    return(paste("'significant' object within 'geneDE' slot must be a",
                 "character with IDs of statiscally significantly",
                 "differentially expressed genes."))
  } else if (!is.data.frame(mirnaTargets(object))) {
    return(paste("'targets' slot must be a data.frame object with miRNAs and",
                 "their relative targets. The end user typically avoids",
                 "manually setting miRNA targets and uses 'getTargets'",
                 "function to retrieve them"))
  } else if (!is.logical(pairedSamples(object))) {
    return(paste("'pairedSamples' must be logical. It should be TRUE if",
                 "miRNA and gene expression data derive from the same samples",
                 "('paired samples') while it should be FALSE if data derive",
                 "from different cohorts of samples"))
  } else if (!is.data.frame(mirnaTargetsIntegration(object))) {
    return(paste("'mirnaTargetsIntegration' must be a data.frame object",
                 "containing miRNA and gene expression data integration.",
                 "The user should use the function 'integrateMirnaTargets()'",
                 "to perform the integration analysis"))
  } else {
    return(TRUE)
  }
})


## -----------
## Constructor
## -----------


#' The constructor function for [MirnaExperiment-class]
#'
#' This is the constructor function that allows to easily create objects of
#' class [`MirnaExperiment`][MirnaExperiment-class]. This function requires as
#' inputs miRNA and gene expression matrices, as well as sample metadata.
#'
#' @details
#' This function requires data to be prepared as described below.
#'
#' ## mirnaExpr and geneExpr
#'
#' `mirnaExpr` and `geneExpr` must be `matrix` objects (or objects coercible
#' to one) that contain miRNA and gene expression values, respectively.
#' Rows must represent the different miRNAs/genes analyzed while columns must
#' represent the different samples in study. For `mirnaExpr`, row names must
#' contain miRNA names according to miRBase nomenclature, whereas for
#' `geneExpr`, row names must contain gene symbols according to hgnc
#' nomenclature. The values contained in these objects can derive from both
#' microarray and RNA-Seq experiments.
#' 
#' For NGS experiments, `mirnaExpr` and `geneExpr` should just be un-normalized
#' count matrices. Instead, for microarray experiments, data should be
#' normalized and log2 transformed, for example with the RMA algorithm.
#'
#' ## samplesMetadata
#'
#' `samplesMetadata` must be a `data.frame` object containing information about
#' samples used for miRNA profiling and for gene expression analysis.
#' Specifically, this data.frame must contain:
#' * A column named `primary`, specifying an identifier for each sample;
#' * A column named `mirnaCol`, containing the column names used for each
#' sample in the `mirnaExpr` object;
#' * A column named `geneCol`, containing the column names used for each
#' sample in the `geneExpr` object;
#' * Other eventual columns that define specific sample metadata, such as
#' disease condition, age, sex and so on...
#' 
#' For unpaired samples, NAs can be used for missing entries in
#' `mirnaCol`/`geneCol`.
#'
#' ## pairedSamples
#'
#' MicroRNA and gene expression measurements may derive from the same subjects
#' (i.e. samples used to generate both miRNA and gene expression data) or from
#' different individuals (i.e. miRNA expression assayed on a group of samples
#' and gene expression retrieved from a different group of samples).
#' `pairedSamples` is a `logical` parameter that defines the relationship
#' between miRNA and gene expression measurements. It must be `TRUE` if data
#' derive from the same individuals, while it must be `FALSE` when data derive
#' from different subjects.
#'
#' @param mirnaExpr A `matrix` object containing microRNA expression levels.
#' Other objects coercible to `matrix` are also accepted (e.g. `data.frame`).
#' This object must be structured as specified in the *details* section
#' @param geneExpr A `matrix` object containing gene expression levels.
#' Other objects coercible to `matrix` are also accepted (e.g. `data.frame`).
#' This object must be structured as specified in the *details* section
#' @param samplesMetadata A `data.frame` object containing information about
#' samples used for microRNA and gene expression profiling. For further
#' information see the *details* section
#' @param pairedSamples Logical, whether miRNA and gene expression levels
#' derive from the same subjects or not. Check the *details* section for
#' additional instructions. Default is `TRUE`
#'
#' @returns
#' A valid [`MirnaExperiment`][MirnaExperiment-class] object containing
#' information about miRNA and gene expression.
#'
#' @examples
#' # load example data
#' data(geneCounts, package = "MIRit")
#' data(mirnaCounts, package = "MIRit")
#' 
#' # create samples metadata
#' meta <- data.frame("primary" = colnames(geneCounts),
#' "mirnaCol" = colnames(mirnaCounts), "geneCol" = colnames(geneCounts),
#' "disease" = c(rep("PTC", 8), rep("NTH", 8)), 
#' "patient" = c(rep(paste("Sample_", seq(8), sep = ""), 2)))
#'
#' # create a 'MirnaExperiment' object
#' obj <- MirnaExperiment(mirnaExpr = mirnaCounts, geneExpr = geneCounts,
#' samplesMetadata = meta, pairedSamples = TRUE)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
MirnaExperiment <- function(
    mirnaExpr,
    geneExpr,
    samplesMetadata,
    pairedSamples = TRUE) {
  
  ## check expression matrices validity
  if (!is.matrix(mirnaExpr) &
      canCoerce(mirnaExpr, "matrix") == FALSE) {
    stop(paste("'mirnaExpr' must be a matrix object or an object coercible",
               "to one. See ?MirnaExperiment for further details on",
               "this object"), call. = FALSE)
  }
  if (!is.matrix(geneExpr) &
      canCoerce(geneExpr, "matrix") == FALSE) {
    stop(paste("'geneExpr' must be a matrix object or an object coercible",
               "to one. See ?MirnaExperiment for further details on",
               "this object"), call. = FALSE)
  }
  if (!is.matrix(mirnaExpr)) {
    mirnaExpr <- as.matrix(mirnaExpr) ## coerce to matrix if needed
  }
  if (!is.matrix(geneExpr)) {
    geneExpr <- as.matrix(geneExpr) ## coerce to matrix if needed
  }
  if (ncol(mirnaExpr) < 2 |
      nrow(mirnaExpr) < 10) {
    stop(paste("'mirnaExpr' must contain expression values deriving from",
               "high-throughput experiments, with samples as columns and",
               "miRNAs as rows. See ?MirnaExperiment for additional details"),
         call. = FALSE)
  }
  if (ncol(geneExpr) < 2 |
      nrow(geneExpr) < 10) {
    stop(paste("'geneExpr' must contain expression values deriving from",
               "high-throughput experiments, with samples as columns and",
               "genes as rows. See ?MirnaExperiment for additional details"),
         call. = FALSE)
  }
  
  ## check metadata provided
  if (!is.data.frame(samplesMetadata) |
      is.null(samplesMetadata$primary) |
      is.null(samplesMetadata$mirnaCol) |
      is.null(samplesMetadata$geneCol) |
      !identical(sort(na.omit(samplesMetadata$mirnaCol)),
                 sort(colnames(mirnaExpr))) |
      !identical(sort(na.omit(samplesMetadata$geneCol)),
                 sort(colnames(geneExpr)))) {
    stop(paste("'samplesMetadata' must be a data.frame object with:\n",
               "\t- one column named 'primary', which contains sample IDs\n",
               "\t- one column named 'mirnaCol', which contains the column",
               "name corresponding to this sample in the expression",
               "table 'mirnaExpr'\n",
               "\t- one column named 'geneCol', which contains the column",
               "name corresponding to this sample in the expression",
               "table 'geneExpr'\n",
               "\t- other columns specifying other information abut samples,",
               "such as age, sex, ecc.\n",
               "For unpaired data, NAs must be used for missing entries",
               "in 'mirnaCol'/'geneCol'"), call. = FALSE)
  }
  
  ## check for valid names in samplesMetadata
  if (any(!colnames(samplesMetadata) %in%
          make.names(colnames(samplesMetadata)))) {
    wrongNames <- which(!colnames(samplesMetadata) %in%
                          make.names(colnames(samplesMetadata)))
    warning(paste("Some variables in the column names of 'samplesMetadata'",
                  "don't have valid R names!", "Therefore,",
                  paste(colnames(samplesMetadata)[wrongNames],
                        collapse = ", "), "will be renamed to:",
                  paste(make.names(colnames(samplesMetadata)[wrongNames]),
                        collapse = ", ")), call. = FALSE)
    colnames(samplesMetadata) <- make.names(colnames(samplesMetadata))
  }
  
  ## check if samples are paired
  if (!is.logical(pairedSamples) |
      length(pairedSamples) != 1) {
    stop(paste("'pairedSamples' must be logical. It should be TRUE if",
               "miRNA and gene expression data derive from the same samples",
               "('paired samples') while it should be FALSE if data derive",
               "from different cohorts of samples"), call. = FALSE)
  }
  if (pairedSamples == TRUE &
      sum(!is.na(samplesMetadata$mirnaCol) &
          !is.na(samplesMetadata$geneCol)) < 3) {
    warning(paste("There are few or no common sample names in",
                  "'samplesMetadata'. Thus, 'pairedSamples' will be set",
                  "to FALSE..."),
            call. = FALSE)
    pairedSamples <- FALSE
  }
  
  ## create a sample map for the MultiAssayExperiment object
  mirnaMap <- samplesMetadata[!is.na(samplesMetadata$mirnaCol), ]
  mirnaMap$geneCol <- NULL
  geneMap <- samplesMetadata[!is.na(samplesMetadata$geneCol), ]
  geneMap$mirnaCol <- NULL
  colnames(mirnaMap)[which(colnames(mirnaMap) == "mirnaCol")] <- "colname"
  colnames(geneMap)[which(colnames(geneMap) == "geneCol")] <- "colname"
  mapList <- list("microRNA" = mirnaMap, "genes" = geneMap)
  sMap <- MultiAssayExperiment::listToMap(mapList)
  
  ## add rownames to metadata table
  rownames(samplesMetadata) <- samplesMetadata$primary
  
  ## create a list with experimental assays
  expList <- list("microRNA" = mirnaExpr, "genes" = geneExpr)
  
  ## create a MultiAssayExperiment object based on user's input
  objMulti <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = expList,
    colData = samplesMetadata,
    sampleMap = sMap)
  
  ## create MirnaExperiment object
  object <- new("MirnaExperiment",
                objMulti,
                pairedSamples = pairedSamples)
  
  ## return the created object
  return(object)
  
}


## ---------
## Accessors
## ---------

#' @rdname mirnaDE
#' @export
setMethod("mirnaDE",
          "MirnaExperiment",
          function(object, onlySignificant, returnObject) {
            if (onlySignificant == TRUE & returnObject == FALSE) {
              object@mirnaDE$data[object@mirnaDE$data$ID %in%
                                    object@mirnaDE$significant, ]
            } else if (onlySignificant == FALSE & returnObject == FALSE){
              object@mirnaDE$data
            } else {
              object@mirnaDE$deObject
            }
          })

#' @rdname geneDE
#' @export
setMethod("geneDE",
          "MirnaExperiment",
          function(object, onlySignificant, returnObject) {
            if (onlySignificant == TRUE & returnObject == FALSE) {
              object@geneDE$data[object@geneDE$data$ID %in%
                                   object@geneDE$significant, ]
            } else if (onlySignificant == FALSE & returnObject == FALSE){
              object@geneDE$data
            } else {
              object@geneDE$deObject
            }
          })

#' @rdname significantMirnas
#' @export
setMethod("significantMirnas", "MirnaExperiment", function(object) {
  object@mirnaDE$significant
})

#' @rdname significantGenes
#' @export
setMethod("significantGenes", "MirnaExperiment", function(object) {
  object@geneDE$significant
})

#' @rdname pairedSamples
#' @export
setMethod("pairedSamples", "MirnaExperiment", function(object) {
  object@pairedSamples
})

#' @rdname mirnaTargets
#' @export
setMethod("mirnaTargets", "MirnaExperiment", function(object) {
  object@targets
})

#' @rdname mirnaTargetsIntegration
#' @export
setMethod("mirnaTargetsIntegration", "MirnaExperiment", function(object) {
  object@mirnaTargetsIntegration
})


## -------
## Setters
## -------

setReplaceMethod("mirnaDE", "MirnaExperiment", function(object, value) {
  object@mirnaDE <- value
  validObject(object)
  object
})

setReplaceMethod("geneDE", "MirnaExperiment", function(object, value) {
  object@geneDE <- value
  validObject(object)
  object
})

setReplaceMethod("pairedSamples", "MirnaExperiment", function(object, value) {
  object@pairedSamples <- value
  validObject(object)
  object
})

setReplaceMethod("mirnaTargets", "MirnaExperiment", function(object, value) {
  object@targets <- value
  validObject(object)
  object
})

setReplaceMethod("mirnaTargetsIntegration",
                 "MirnaExperiment",
                 function(object, value) {
                   object@mirnaTargetsIntegration <- value
                   validObject(object)
                   object
                 })



## ==========================================================================
## MirnaEnrichment Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The `MirnaEnrichment` class
#'
#' This class extends and adapts the
#' [`enrichResult-class`][DOSE::enrichResult-class] in order to make it
#' suitable for handling miRNA enrichment results.
#'
#' @slot result A `data.frame` object holding the output of enrichment analysis
#' @slot pvalueCutoff A `numeric` value defining the threshold used for
#' statistical significance in the enrichment analysis (e.g. `0.05`)
#' @slot pAdjustMethod A `character` indicating the method used to correct
#' p-values for multiple testing (e.g. `fdr`)
#' @slot qvalueCutoff A `numeric` value defining the q-value threshold used
#' in the enrichment analysis (e.g. `0.2`)
#' @slot organism The name of the organism under consideration (e.g.
#' `Homo sapiens`)
#' @slot ontology The name of the miEAA 2.0 category used to perform the
#' enrichment analysis (e.g. `miRWalk_GO_mature`)
#' @slot gene A `character` vector containing the list of miRNAs used for the
#' enrichment
#' @slot keytype The type of miRNA IDs used. For example `miRBase v22`
#' @slot universe The background universe of miRNAs used for the
#' over-representation analysis (ORA). Typically, this is equal to the complete
#' list of miRNAs assayed
#' @slot gene2Symbol Mapping of genes to symbols, if needed
#' @slot geneSets Gene sets
#' @slot readable Logical flag of gene ID in symbol or not
#' @slot termsim Similarity between terms
#' @slot method Method for calculating the similarity between nodes
#' @slot dr Dimension reduction result
#' 
#' @param object An object of class [`MirnaEnrichment`][MirnaEnrichment-class]
#' @param x An object of class [`MirnaEnrichment`][MirnaEnrichment-class]
#' containing enrichment results of multiple miEAA 2.0 categories
#' @param i A valid miEAA 2.0 category name
#' @param j Missing
#' @param drop Missing
#' @param ... Additional arguments not used
#'
#' @references
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name MirnaEnrichment-class
#' @docType class
#' @export
#' @import methods
#' @importClassesFrom DOSE enrichResult
setClass("MirnaEnrichment",
         contains="enrichResult")


## ---------
## Accessors
## ---------

#' @rdname enrichmentResults
#' @export
setMethod("enrichmentResults", "MirnaEnrichment", function(object) {
  object@result
})

setMethod("enrichmentDatabase", "MirnaEnrichment", function(object) {
  object@ontology
})

setMethod("mirnaIdEnrichment", "MirnaEnrichment", function(object) {
  object@gene
})


## -------
## Setters
## -------

setReplaceMethod("enrichmentResults", "MirnaEnrichment", function(object, value) {
  object@result <- value
  validObject(object)
  object
})

setReplaceMethod("enrichmentDatabase", "MirnaEnrichment", function(object, value) {
  object@ontology <- value
  validObject(object)
  object
})



## ==========================================================================
## MirnaGsea Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The `MirnaGsea` class
#'
#' This class extends and adapts the
#' [`gseaResult-class`][DOSE::gseaResult-class] in order to make it
#' suitable for handling miRNA gene set enrichment analysis (GSEA) results.
#'
#' @slot result A `data.frame` object holding the output of enrichment analysis
#' @slot pvalueCutoff A `numeric` value defining the threshold used for
#' statistical significance in the enrichment analysis (e.g. `0.05`)
#' @slot pAdjustMethod A `character` indicating the method used to correct
#' p-values for multiple testing (e.g. `fdr`)
#' @slot organism The name of the organism under consideration (e.g.
#' `Homo sapiens`)
#' @slot ontology The name of the miEAA 2.0 category used to perform the
#' enrichment analysis (e.g. `miRWalk_GO_mature`)
#' @slot gene A `character` vector containing the list of miRNAs used for the
#' enrichment
#' @slot lfc A `numeric` vector containing the list of log2 fold changes
#' relative to miRNAs listed in the `gene` slot
#' @slot keytype The type of miRNA IDs used. For example `miRBase v22`
#' @slot gene2Symbol Mapping of genes to symbols, if needed
#' @slot geneSets Gene sets
#' @slot setType Set type
#' @slot geneList Order rank gene list
#' @slot permScores Permutation scores
#' @slot params Parameters
#' @slot readable Logical flag of gene ID in symbol or not
#' @slot dr Dimension reduction result
#' 
#' @param object An object of class [`MirnaGsea`][MirnaGsea-class]
#' @param x An object of class [`MirnaGsea`][MirnaGsea-class] containing
#' enrichment results of multiple miEAA 2.0 categories
#' @param i A valid miEAA 2.0 category name
#' @param j Missing
#' @param drop Missing
#' @param ... Additional arguments not used
#'
#' @references
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name MirnaGsea-class
#' @docType class
#' @export
#' @import methods
#' @importClassesFrom DOSE gseaResult
setClass("MirnaGsea",
         contains="gseaResult",
         slots = representation(
           pvalueCutoff = "numeric",
           pAdjustMethod = "character",
           ontology = "character",
           gene = "character",
           lfc = "numeric"
         ))


## ---------
## Accessors
## ---------

#' @rdname enrichmentResults
#' @export
setMethod("enrichmentResults", "MirnaGsea", function(object) {
  object@result
})

setMethod("enrichmentDatabase", "MirnaGsea", function(object) {
  object@ontology
})

setMethod("lfcEnrichment", "MirnaGsea", function(object) {
  object@lfc
})

setMethod("mirnaIdEnrichment", "MirnaGsea", function(object) {
  object@gene
})


## -------
## Setters
## -------

setReplaceMethod("enrichmentResults", "MirnaGsea", function(object, value) {
  object@result <- value
  validObject(object)
  object
})

setReplaceMethod("enrichmentDatabase", "MirnaGsea", function(object, value) {
  object@ontology <- value
  validObject(object)
  object
})

