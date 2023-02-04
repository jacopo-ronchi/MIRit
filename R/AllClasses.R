## ==========================================================================
## MirnaExperiment Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The 'MirnaExperiment' class
#'
#' @description
#' This class extends the [`MultiAssayExperiment::MultiAssayExperiment-class`]
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
#' @slot ExperimentList A [`MultiAssayExperiment::ExperimentList`] class object
#' for each assay dataset
#' @slot colData A `DataFrame` of all clinical/specimen data available across
#' experiments
#' @slot sampleMap A `DataFrame` of translatable identifiers of samples and
#' participants
#' @slot metadata Additional data describing the object
#' @slot drops A metadata `list` of dropped information
#' @slot mirnaDE A `data.frame` object containing miRNA differential expression
#' results
#' @slot geneDE A `data.frame` object containing gene differential expression
#' results
#' @slot significantMirnas A `character` vector containing the names of
#' significantly differentially expressed miRNAs
#' @slot significantGenes A `character` vector containing the names of
#' significantly differentially expressed genes
#' @slot pairedSamples A `logical` parameter that specifies whether miRNA and
#' gene expression measurements derive from the same individuals (`TRUE`) or
#' from different subjects (`FALSE`)
#' @slot targets A `data.frame` object containing miRNA-target pairs. This
#' slot is commonly populated by the [getTargets()] function
#' @slot mirnaTargetsIntegration A `data.frame` object containing the results
#' of the integration analysis between miRNA and gene expression values. This
#' slot is commonly populated by the [integrateMirnaTargets()] function
#'
#' @inheritSection
#' MultiAssayExperiment::MultiAssayExperiment-class ExperimentList
#'
#' @inheritSection MultiAssayExperiment::MultiAssayExperiment-class colData
#'
#' @inheritSection MultiAssayExperiment::MultiAssayExperiment-class sampleMap
#'
#' # mirnaDE and geneDE
#'
#' `mirnaDE` consists of a `data.frame` with five columns:
#' * `ID`: indicates the name of the miRNA;
#' * `logFC`: indicates the fold change of each miRNA in logarithmic scale;
#' * `AveExpr`: represents the average expression of each miRNA;
#' * `P.Value`: indicates the p-value resulting from differential expression
#' analysis;
#' * `FDR`: contains the p-values adjusted for multiple testing.
#'
#' `geneDE` consists of a `data.frame` with five columns:
#' * `ID`: indicates the name of the gene;
#' * `logFC`: indicates the fold change of each gene in logarithmic scale;
#' * `AveExpr`: represents the average expression of each gene;
#' * `P.Value`: indicates the p-value resulting from differential expression
#' analysis;
#' * `FDR`: contains the p-values adjusted for multiple testing.
#'
#' # significantMirnas and significantGenes
#'
#' `significantMirnas` and `significantGenes` are `character` vectors of
#' miRNA/gene IDs that are considered as statistically differentially expressed.
#' The IDs contained in these vectors must be present in the `ID` columns of
#' `mirnaDE` and `geneDE`.
#'
#' # pairedSamples
#'
#' As already mentioned, `pairedSamples` must be `TRUE` when miRNA and gene
#' expression derive from the same subjects, while it must `FALSE` if this is
#' not the case.
#'
#' # targets
#'
#' `targets` is a `data.frame` with just two columns:
#' * `mature_mirna_id`, which contains miRNA names; and
#' * `target_symbol`, which indicates the target gene for the corresponding
#' miRNA.
#'
#' # mirnaTargetsIntegration
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
#' * `target`: the correlated target gene;
#' * `microRNA.Direction`: the fold change direction of the DE-miRNA;
#' * `Correlation.Type`: specifies the correlation direction, either `positive`
#' or `negative`;
#' * `Corr.Coefficient`: the value of the correlation coefficient used;
#' * `Corr.P.Value`: the p-value resulting from the correlation analysis;
#' * `Corr.Adjusted.P.Val`: contains the correlation p-values corrected for
#' multiple testing.
#'
#' @note
#' To create a [`MirnaExperiment`][MirnaExperiment-class] object, you can use
#' the [MirnaExperiment()] constructor function, which allows to easily build
#' and verify a valid object starting from the results of differential
#' expression analysis.
#'
#' @name MirnaExperiment-class
#'
#' @references
#' Marcel Ramos et al. Software For The Integration Of Multiomics Experiments
#' In Bioconductor. Cancer Research, 2017 November 1; 77(21); e39-42. DOI:
#' [10.1158/0008-5472.CAN-17-0344]
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @seealso
#' See [MultiAssayExperiment::MultiAssayExperiment-class] for additional
#' information.
#'
#' @docType class
#' @export
#' @import methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setClass("MirnaExperiment",
         contains="MultiAssayExperiment",
         slots=representation(
           mirnaDE="data.frame",
           significantMirnas="character",
           geneDE="data.frame",
           significantGenes="character",
           pairedSamples="logical",
           targets="data.frame",
           mirnaTargetsIntegration="data.frame"))


## --------
## Validity
## --------

setValidity("MirnaExperiment", function(object) {

  if (!is.data.frame(mirnaDE(object))) {
    return(paste("'mirnaDE' slot must be a data.frame object with miRNA",
                 "differential expression results, such as the output of",
                 "'topTable' in limma"))
  } else if (!is.data.frame(geneDE(object))) {
    return(paste("'geneDE' slot must be a data.frame object with gene",
                 "differential expression results, such as the output of",
                 "'topTable' in limma"))
  } else if (!identical(colnames(mirnaDE(object)), c("ID", "logFC", "AveExpr", "P.Value", "FDR")) |
             !identical(colnames(geneDE(object)), c("ID", "logFC", "AveExpr", "P.Value", "FDR"))) {
    return(paste("'mirnaDE' and 'geneDE' slots must be dataframes with",
                 "column names: 'ID', 'logFC', 'AveExpr', 'P.Value', 'FDR'."))
  } else if (!is.character(significantMirnas(object))  |
             !all(significantMirnas(object) %in% mirnaDE(object,
                                                         onlySignificant = FALSE)$ID)) {
    return(paste("'significantMirnas' must be a character with IDs of",
                 "statiscally significantly differentially expressed miRNAs."))
  } else if (!is.character(significantGenes(object))  |
             !all(significantGenes(object) %in% geneDE(object,
                                                       onlySignificant = FALSE)$ID)) {
    return(paste("'significantGenes' must be a character with IDs of",
                 "statiscally significantly differentially expressed genes"))
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
#' inputs miRNA and gene expression matrices as well as the results of miRNA
#' and gene differential expression analysis.
#'
#' @details
#' This function requires data to be prepared as described below.
#'
#' ## mirnaExpr and geneExpr
#'
#' `mirnaExpr` and `geneExpr` must be `matrix` objects (or objects coercible
#' to one) that contain miRNA and gene expression values, respectively.
#' Rows must represent the different miRNAs/genes assayed while columns must
#' represent the different samples analyzed. For `mirnaExpr`, row names must
#' contain miRNA names according to miRBase nomenclature, whereas for
#' `geneExpr`, row names must contain gene symbols according to hgnc
#' nomenclature. The values contained in these objects can derive from both
#' microarray and RNA-Seq experiments. For microarray experiments, data should
#' be normalized and log2 transformed, for example with the RMA algorithm,
#' while for NGS experiments, data should consist of normalized counts/CPM.
#' For instance, normalized values for NGS experiments can be extracted from a
#' [DESeq2] object with `counts(obj, normalized = TRUE)`, or from an [edgeR]
#' object through `cpm(obj, normalized.lib.sizes = TRUE)`.
#'
#' ## mirnaMetadata and geneMetadata
#'
#' `mirnaMetadata` and `geneMetadata` must be two `data.frame` objects
#' containing information about samples used for miRNA profiling and for gene
#' expression analysis, respectively. Specifically they must contain:
#' * A column named `primary`, specifying an identifier for each sample;
#' * A column named `colName`, containing the column names used for each sample
#' in the `mirnaExpr`/`geneExpr` object;
#' * Other eventual columns that define specific sample metadata, such as age,
#' sex and so on...
#'
#' ## mirnaDE and geneDE
#'
#' `mirnaDE` and `geneDE` are two objects of class `data.frame` containing
#' the results of miRNA and gene differential expression analysis respectively.
#' These tables should contain the differential expression results for all
#' miRNAs/genes analyzed, not just for statistically significant species.
#' This function directly supports the output of the most common functions
#' used by different packages for differential expression, such as
#' [limma::topTable()], [DESeq2::results()] and [edgeR::topTags()].
#' However, other `data.frame` objects can be used, as long as they have:
#' * One column containing miRNA/gene names (according to miRBase/hgnc
#' nomenclature). Accepted column names are: `ID`, `Symbol`, `Gene_Symbol`,
#' `Mirna`, `mir`, `Gene`, `gene.symbol`, `Gene.symbol`;
#' * One column with log2 fold changes. Accepted column names are: `logFC`,
#' `log2FoldChange`, `FC`, `lFC`;
#' * One column with average expression. Accepted column names are: `AveExpr`,
#' `baseMean`, `logCPM`;
#' * One column with the p-values resulting from the differential expression
#' analysis. Accepted column names are: `P.Value`, `pvalue`, `PValue`,
#' `Pvalue`;
#' * One column containing p-values adjusted for multiple testing. Accepted
#' column names are: `adj.P.Val`, `padj`, `FDR`, `fdr`, `adj`, `adj.p`, `adjp`.
#'
#' ## significantMirnas and significantGenes
#'
#' `significantMirnas` and `significantGenes` are two `character` vectors that
#' specifies the IDs of miRNAs and genes considered to be significantly
#' differentially expressed. The miRNA IDs contained in `significantMirnas`
#' must be present in the `ID` column of `mirnaDE`, in the same way as gene
#' symbols in `significantGenes` must be present in the `ID` column
#' of `geneDE`.
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
#' from different subjects. It is important to note that, for paired samples,
#' the sample IDs contained in the `primary` column of `mirnaMetadata` should
#' match those contained in the `primary` column of `geneMetadata`.
#'
#' @param mirnaExpr A `matrix` object containing microRNA expression levels.
#' Other objects coercible to `matrix` are also accepted (e.g. `data.frame`).
#' This object must be structured as specified in the *details* section
#' @param geneExpr A `matrix` object containing gene expression levels.
#' Other objects coercible to `matrix` are also accepted (e.g. `data.frame`).
#' This object must be structured as specified in the *details* section
#' @param mirnaMetadata A `data.frame` object containing information about
#' samples used for microRNA profiling. For further information see the
#' *details* section
#' @param geneMetadata A `data.frame` object containing information about
#' samples used for gene expression analysis. For further information see the
#' *details* section
#' @param mirnaDE A `data.frame` containing the output of miRNA differential
#' expression analysis. Check the *details* section to see the required format
#' @param geneDE A `data.frame` containing the output of gene differential
#' expression analysis. Check the *details* section to see the required format
#' @param significantMirnas A `character` vector containing the IDs of
#' statistically differentially expressed miRNAs. See the *details* section for
#' further information
#' @param significantGenes A `character` vector containing the IDs of
#' statistically differentially expressed genes. See the *details* section for
#' further information
#' @param Logical, wheteher miRNA and gene expression levels derive from the
#' same subjects or not. Check the *details* section for additional
#' instructions. Default is `TRUE`
#'
#' @returns
#' A valid [`MirnaExperiment`][MirnaExperiment-class] object containing
#' information about miRNA and gene expression.
#'
#' @examples
#' # limma output for miRNAs
#' de_m <- limma::topTable(fit_m, number = Inf)
#' sig_m <- de_m[1, which(de_m$adj.P.Val < 0.05)]
#'
#' # limma output for genes
#' de_g <- limma::topTable(fit_g, number = Inf)
#' sig_g <- de_g[1, which(de_g$adj.P.Val < 0.05)]
#'
#' # create a 'MirnaExperiment' object after DE analysis with limma
#' obj <- MirnaExperiment(mirnaExpr = exp_m, geneExpr = exp_b,
#' mirnaMetadata = info_m, geneMetadata = info_g, mirnaDE = de_m,
#' geneDE = de_g, significantMirnas = sig_m, significantGenes = sig_g,
#' pairedSamples = TRUE)
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
MirnaExperiment <- function(
    mirnaExpr,
    geneExpr,
    mirnaMetadata,
    geneMetadata,
    mirnaDE,
    geneDE,
    significantMirnas,
    significantGenes,
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

  ## check differential expression inputs
  if (!is.data.frame(mirnaMetadata) |
      is.null(mirnaMetadata$primary) |
      is.null(mirnaMetadata$colName) |
      !identical(sort(mirnaMetadata$colName), sort(colnames(mirnaExpr)))) {
    stop(paste("'mirnaMetadata' must be a data.frame object with:\n",
               "\t- one column named 'primary', which contains the sample ID\n",
               "\t- one column named 'colName', which contains the column",
               "name corresponding to this sample in the expression",
               "table 'mirnaExpr'\n",
               "\t- other columns specifying other information abut samples,",
               "such as age, sex, ecc."), call. = FALSE)
  }
  if (!is.data.frame(geneMetadata) |
      is.null(geneMetadata$primary) |
      is.null(geneMetadata$colName) |
      !identical(sort(geneMetadata$colName), sort(colnames(geneExpr)))) {
    stop(paste("'geneMetadata' must be a data.frame object with:\n",
               "\t- one column named 'primary', which contains the sample ID\n",
               "\t- one column named 'colName', which contains the column",
               "name corresponding to this sample in the expression",
               "table 'geneExpr'\n",
               "\t- other columns specifying other information abut samples,",
               "such as age, sex, ecc."), call. = FALSE)
  }
  if (!is.data.frame(mirnaDE)) {
    stop(paste("'mirnaDE' slot must be a data.frame object with miRNA",
               "differential expression results, such as the output of",
               "'topTable' in limma"), call. = FALSE)
  }
  if (!is.data.frame(geneDE)) {
    stop(paste("'geneDE' slot must be a data.frame object with gene",
               "differential expression results, such as the output of",
               "'topTable' in limma"), call. = FALSE)
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
      length(intersect(mirnaMetadata$primary, geneMetadata$primary)) < 3) {
    warning(paste("There are few or no common sample names in the 'primary'",
                  "columns of 'mirnaMetadata' and 'geneMetadata'.",
                  "Thus, pairedSamples will be set to FALSE..."),
            call. = FALSE)
    pairedSamples <- FALSE
  }

  ## check and set accepted column names in DE data.frames
  mirnaDE <- identifyColNames(mirnaDE, "miRNA")
  geneDE <- identifyColNames(geneDE, "gene")

  ## check that miRNA and gene names are the same across data
  if (!identical(sort(rownames(mirnaExpr)), sort(mirnaDE$ID))) {
    stop(paste("Row names of 'mirnaExpr' data.frame must match the miRNA",
               "identifiers present in 'mirnaDE'"), call. = FALSE)
  }
  if (!identical(sort(rownames(geneExpr)), sort(geneDE$ID))) {
    stop(paste("Row names of 'geneExpr' data.frame must match the gene",
               "identifiers present in 'geneDE'"), call. = FALSE)
  }

  ## check that significant miRNAs/genes are a subset of all miRNA/genes tested
  if (!is.character(significantMirnas) |
      !all(significantMirnas %in% mirnaDE$ID)) {
    stop(paste("'significantMirnas' must be a character with IDs of",
               "statiscally significantly differentially expressed miRNAs.",
               "They must match the identifiers present in 'mirnaDE'."),
         call. = FALSE)
  }
  if (!is.character(significantGenes) |
      !all(significantGenes %in% geneDE$ID)) {
    stop(paste("'significantGenes' must be a character with IDs of",
               "statiscally significantly differentially expressed genes.",
               "They must match the identifiers present in 'geneDE'."),
         call. = FALSE)
  }

  ## create a MultiAssayExperiment object based on user's input
  mirnaMap <- data.frame("primary"=mirnaMetadata$primary,
                         "colname"=mirnaMetadata$colName)
  genesMap <- data.frame("primary"=geneMetadata$primary,
                         "colname"=geneMetadata$colName)
  mapList <- list("microRNA"=mirnaMap, "genes"=genesMap)
  sMap <- listToMap(mapList)
  colnames(mirnaMetadata)[which(colnames(mirnaMetadata)
                                == "colName")] <- "mirnaColName"
  colnames(geneMetadata)[which(colnames(geneMetadata)
                                == "colName")] <- "geneColName"
  samplesMetadata <- suppressMessages(dplyr::full_join(mirnaMetadata,
                                                       geneMetadata,
                                                       by = "primary"))
  rownames(samplesMetadata) <- samplesMetadata$primary
  expList <- list("microRNA"=mirnaExpr, "genes"=geneExpr)
  objMulti <- MultiAssayExperiment(experiments=expList,
                                   colData=samplesMetadata,
                                   sampleMap=sMap)

  ## create MirnaExperiment object
  object <- new("MirnaExperiment",
                objMulti,
                mirnaDE=mirnaDE,
                geneDE=geneDE,
                significantMirnas=significantMirnas,
                significantGenes=significantGenes,
                pairedSamples=pairedSamples)

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
          function(object, onlySignificant) {
            if (onlySignificant == TRUE) {
              object@mirnaDE[object@mirnaDE$ID %in% object@significantMirnas, ]
            } else {
              object@mirnaDE
            }
          })

#' @rdname geneDE
#' @export
setMethod("geneDE",
          "MirnaExperiment",
          function(object, onlySignificant) {
            if (onlySignificant == TRUE) {
              object@geneDE[object@geneDE$ID %in% object@significantGenes, ]
            } else {
              object@geneDE
            }
          })

#' @rdname significantMirnas
#' @export
setMethod("significantMirnas", "MirnaExperiment", function(object) {
  object@significantMirnas
})

#' @rdname significantGenes
#' @export
setMethod("significantGenes", "MirnaExperiment", function(object) {
  object@significantGenes
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

setReplaceMethod("significantMirnas", "MirnaExperiment", function(object, value) {
  object@significantMirnas <- value
  validObject(object)
  object
})

setReplaceMethod("significantGenes", "MirnaExperiment", function(object, value) {
  object@significantGenes <- value
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
#' This function extends and adapts the [`DOSE::enrichResult-class`] in order
#' to make it suitable for handling miRNA enrichment results.
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
#' @slot termism Similarity between terms
#' @slot method Method for calculating the similarity between nodes
#' @slot dr Dimension reduction result
#'
#' @references
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
#' @import methods
#' @importClassesFrom clusterProfiler enrichResult
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
#' This function extends and adapts the [`DOSE::gseaResult-class`] in order
#' to make it suitable for handling miRNA gene set enrichment analysis
#' (GSEA) results.
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
#' @references
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
#' @import methods
#' @importClassesFrom clusterProfiler gseaResult
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

