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
#' @slot integration A `list` object containing the results
#' of the integration analysis between miRNA and gene expression values. This
#' slot is commonly populated by the [mirnaIntegration()] function
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
#' * `group`, which is the column name of the variable (in colData) used for
#' differential expression analysis;
#' * `contrast`, which represents the groups that are compared during
#' differential expression analysis (e.g. 'disease-healthy');
#' * `design`, which outlines the R `formula` used for fitting the model. It
#' includes the variable of interest (`group`) together with eventual
#' covariates (e.g. '~ 0 + disease + sex');
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
#' `targets` is a `data.frame` with miRNA-target interactions, as retrieved by
#' [getTargets()] function.
#'
#' @section integration:
#'
#' Lastly, `integration` slot contains a `list` object that stores
#' the results and the options used for performing the integrative miRNA-gene
#' analysis. In particular, `integration` contains:
#' * `data`, which is a `data.frame` object with the results of the integrative
#' analysis;
#' * `method`, which specifies the procedure used to perform the integrative
#' analysis;
#' * `pCutoff`, which indicates the p-value cutoff used for the analysis;
#' * `pAdjustment`, the approach used for multiple testing correction.
#'
#' Moreover, `data` differs on the basis of the integration strategy used. For
#' the one-sided association test integration, and for integration based on
#' rotation gene set tests, this `data.frame` has seven columns:
#' * `microRNA`: the miRNA ID;
#' * `mirna.direction`: the fold change direction of the DE-miRNA (`Up` or
#' `Down`);
#' * `gene.direction`: the fold change direction of target genes (`Up` or
#' `Down`);
#' * `DE`: represents the number of differentially expressed targets;
#' * `targets`: represents the total number of targets for this miRNA;
#' * `P.Val`: indicates the resulting p-value;
#' * `adj.P.Val`: contains test p-values corrected for multiple testing;
#' * `DE.targets`: contains the list of differentially expressed targets whose
#' expression is negatively associated with miRNA expression.
#'
#' Instead, when a correlation analysis is performed, `data` has six columns:
#' * `microRNA`: the miRNA ID;
#' * `Target`: the correlated target gene;
#' * `microRNA.Direction`: the fold change direction of the DE-miRNA;
#' * `Corr.Coeff`: the value of the correlation coefficient used;
#' * `Corr.P.Value`: the p-value resulting from the correlation analysis;
#' * `Corr.Adjusted.P.Val`: contains the correlation p-values corrected for
#' multiple testing.
#'
#' To access the results of the integrative analysis, the `data` slot can be
#' accessed through the [integration()] function.
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
#' @import MultiAssayExperiment
setClass("MirnaExperiment",
         contains = "MultiAssayExperiment",
         slots = representation(
             mirnaDE = "list",
             geneDE = "list",
             pairedSamples = "logical",
             targets = "data.frame",
             integration = "list"
         )
)


## -----------------
## Initialize method
## -----------------


setMethod(
    "initialize",
    signature(.Object = "MirnaExperiment"),
    function(.Object, ...) {
        deList <- list(
            data = data.frame(),
            significant = character(),
            method = character(),
            group = character(),
            contrast = character(),
            design = NULL,
            pCutoff = numeric(),
            pAdjustment = character(),
            logFC = numeric(),
            deObject = NULL
        )
        
        intList <- list(
            data = data.frame(),
            method = character(),
            pCutoff = numeric(),
            pAdjustment = character()
        )
        
        .Object <- callNextMethod(.Object,
                                  ...,
                                  mirnaDE = deList,
                                  geneDE = deList,
                                  integration = intList
        )
        .Object
    }
)


## --------
## Validity
## --------

setValidity("MirnaExperiment", function(object) {
    if (!is.list(object@mirnaDE)) {
        return(paste(
            "'mirnaDE' slot must be a list object with miRNA",
            "differential expression results. Please see",
            "?MirnaExperiment-class"
        ))
    } else if (!is.list(object@geneDE)) {
        return(paste(
            "'geneDE' slot must be a list object with gene",
            "differential expression results. Please see",
            "?MirnaExperiment-class"
        ))
    } else if (!identical(
        sort(names(object@mirnaDE)),
        sort(c(
            "data", "significant", "method", "group",
            "contrast", "design", "pCutoff",
            "pAdjustment", "logFC", "deObject"
        ))
    ) |
    !identical(
        sort(names(object@geneDE)),
        sort(c(
            "data", "significant", "method", "group",
            "contrast", "design", "pCutoff",
            "pAdjustment", "logFC", "deObject"
        ))
    )) {
        return(paste(
            "'mirnaDE' and 'geneDE' slots must be list objects",
            "containing: 'data', 'significant', 'method', 'group',",
            "'contrast', 'design', 'pCutoff', 'pAdjustment', 'logFC',",
            "and 'deObject'. Please see ?MirnaExperiment-class"
        ))
    } else if (!is.data.frame(mirnaDE(object))) {
        return(paste(
            "'data' within mirnaDE slot must be a data.frame object with",
            "miRNA differential expression results. Please see",
            "?MirnaExperiment-class"
        ))
    } else if (!is.data.frame(geneDE(object))) {
        return(paste(
            "'data' within geneDE slot must be a data.frame object with",
            "gene differential expression results. Please see",
            "?MirnaExperiment-class"
        ))
    } else if (!is.character(significantMirnas(object)) |
               !all(significantMirnas(object) %in%
                    mirnaDE(object, onlySignificant = FALSE)$ID)) {
        return(paste(
            "'significant' object within 'mirnaDE' slot must be a",
            "character with IDs of statiscally significantly",
            "differentially expressed miRNAs."
        ))
    } else if (!is.character(significantGenes(object)) |
               !all(significantGenes(object)
                    %in% geneDE(object, onlySignificant = FALSE)$ID)) {
        return(paste(
            "'significant' object within 'geneDE' slot must be a",
            "character with IDs of statiscally significantly",
            "differentially expressed genes."
        ))
    } else if (!is.data.frame(mirnaTargets(object))) {
        return(paste(
            "'targets' slot must be a data.frame object with miRNAs and",
            "their relative targets. The end user typically avoids",
            "manually setting miRNA targets and uses 'getTargets'",
            "function to retrieve them."
        ))
    } else if (!is.logical(pairedSamples(object))) {
        return(paste(
            "'pairedSamples' must be logical. It should be TRUE if",
            "miRNA and gene expression data derive from the same samples",
            "('paired samples') while it should be FALSE if data derive",
            "from different cohorts of samples"
        ))
    } else if (!is.list(object@integration)) {
        return(paste(
            "'integration' slot must be a list object with",
            "integrative analysis results. Please see",
            "?MirnaExperiment-class"
        ))
    } else if (!is.data.frame(integration(object))) {
        return(paste(
            "'data' within integration slot must be a data.frame",
            "containing miRNA and gene expression data integration.",
            "The user should use the function 'mirnaIntegration()'",
            "to perform the integration analysis"
        ))
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
#' meta <- data.frame(
#'     "primary" = colnames(geneCounts),
#'     "mirnaCol" = colnames(mirnaCounts), "geneCol" = colnames(geneCounts),
#'     "disease" = c(rep("PTC", 8), rep("NTH", 8)),
#'     "patient" = c(rep(paste("Sample_", seq(8), sep = ""), 2))
#' )
#'
#' # create a 'MirnaExperiment' object
#' obj <- MirnaExperiment(
#'     mirnaExpr = mirnaCounts, geneExpr = geneCounts,
#'     samplesMetadata = meta, pairedSamples = TRUE
#' )
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
MirnaExperiment <- function(mirnaExpr,
                            geneExpr,
                            samplesMetadata,
                            pairedSamples = TRUE) {
    ## check expression matrices validity
    if (!is.matrix(mirnaExpr) &
        canCoerce(mirnaExpr, "matrix") == FALSE) {
        stop("'mirnaExpr' must be a matrix object or an object coercible ",
             "to one. See ?MirnaExperiment for further details on ",
             "this object",
             call. = FALSE
        )
    }
    if (!is.matrix(geneExpr) &
        canCoerce(geneExpr, "matrix") == FALSE) {
        stop("'geneExpr' must be a matrix object or an object coercible ",
             "to one. See ?MirnaExperiment for further details on ",
             "this object",
             call. = FALSE
        )
    }
    if (!is.matrix(mirnaExpr)) {
        mirnaExpr <- as.matrix(mirnaExpr) ## coerce to matrix if needed
    }
    if (!is.matrix(geneExpr)) {
        geneExpr <- as.matrix(geneExpr) ## coerce to matrix if needed
    }
    if (ncol(mirnaExpr) < 2 |
        nrow(mirnaExpr) < 10) {
        stop("'mirnaExpr' must contain expression values deriving from ",
             "high-throughput experiments, with samples as columns and ",
             "miRNAs as rows. See ?MirnaExperiment for additional details",
             call. = FALSE
        )
    }
    if (ncol(geneExpr) < 2 |
        nrow(geneExpr) < 10) {
        stop("'geneExpr' must contain expression values deriving from ",
             "high-throughput experiments, with samples as columns and ",
             "genes as rows. See ?MirnaExperiment for additional details",
             call. = FALSE
        )
    }
    
    ## check metadata provided
    if (!is.data.frame(samplesMetadata) |
        is.null(samplesMetadata$primary) |
        is.null(samplesMetadata$mirnaCol) |
        is.null(samplesMetadata$geneCol) |
        !identical(
            sort(na.omit(samplesMetadata$mirnaCol)),
            sort(colnames(mirnaExpr))
        ) |
        !identical(
            sort(na.omit(samplesMetadata$geneCol)),
            sort(colnames(geneExpr))
        )) {
        stop("'samplesMetadata' must be a data.frame object with:\n",
             "\t- one column named 'primary', which contains sample IDs\n",
             "\t- one column named 'mirnaCol', which contains the column ",
             "name corresponding to this sample in the expression ",
             "table 'mirnaExpr'\n",
             "\t- one column named 'geneCol', which contains the column ",
             "name corresponding to this sample in the expression ",
             "table 'geneExpr'\n",
             "\t- other columns specifying other information abut samples, ",
             "such as age, sex, ecc.\n",
             "For unpaired data, NAs must be used for missing entries ",
             "in 'mirnaCol'/'geneCol'",
             call. = FALSE
        )
    }
    
    ## check for valid names in samplesMetadata
    if (any(!colnames(samplesMetadata) %in%
            make.names(colnames(samplesMetadata)))) {
        wrongNames <- which(!colnames(samplesMetadata) %in%
                                make.names(colnames(samplesMetadata)))
        warning("Some variables in the column names of 'samplesMetadata' ",
                "don't have valid R names! ", "Therefore, ",
                paste(colnames(samplesMetadata)[wrongNames],
                      collapse = ", "
                ), " will be renamed to: ",
                paste(make.names(colnames(samplesMetadata)[wrongNames]),
                      collapse = ", "
                ),
                call. = FALSE
        )
        colnames(samplesMetadata) <- make.names(colnames(samplesMetadata))
    }
    
    ## check if samples are paired
    if (!is.logical(pairedSamples) |
        length(pairedSamples) != 1) {
        stop("'pairedSamples' must be logical. It should be TRUE if ",
             "miRNA and gene expression data derive from the same samples ",
             "('paired samples') while it should be FALSE if data derive ",
             "from different cohorts of samples",
             call. = FALSE
        )
    }
    if (pairedSamples == TRUE &
        sum(!is.na(samplesMetadata$mirnaCol) &
            !is.na(samplesMetadata$geneCol)) < 3) {
        warning("There are few or no common sample names in ",
                "'samplesMetadata'. Thus, 'pairedSamples' will be set ",
                "to FALSE...",
                call. = FALSE
        )
        pairedSamples <- FALSE
    }
    
    ## check that primary column doesn't have repeated elements (unpaired data)
    if (pairedSamples == FALSE &
        any(duplicated(samplesMetadata$primary)) == TRUE) {
        stop("For unpaired data, 'primary' column in 'samplesMetadata' ",
             "can't contain duplicated elements. See ?MirnaExperiment",
             call. = FALSE
        )
    }
    
    ## create a sample map for the MultiAssayExperiment object
    mirnaMap <- samplesMetadata[!is.na(samplesMetadata$mirnaCol), ]
    mirnaMap$geneCol <- NULL
    geneMap <- samplesMetadata[!is.na(samplesMetadata$geneCol), ]
    geneMap$mirnaCol <- NULL
    colnames(mirnaMap)[which(colnames(mirnaMap) == "mirnaCol")] <- "colname"
    colnames(geneMap)[which(colnames(geneMap) == "geneCol")] <- "colname"
    geneMap <- geneMap[, colnames(mirnaMap)]
    mapList <- list("microRNA" = mirnaMap, "genes" = geneMap)
    sMap <- listToMap(mapList)
    
    ## add rownames to metadata table
    rownames(samplesMetadata) <- samplesMetadata$primary
    
    ## create a list with experimental assays
    expList <- list("microRNA" = mirnaExpr, "genes" = geneExpr)
    
    ## create a MultiAssayExperiment object based on user's input
    objMulti <- MultiAssayExperiment(
        experiments = expList,
        colData = samplesMetadata,
        sampleMap = sMap
    )
    
    ## create MirnaExperiment object
    object <- new("MirnaExperiment",
                  objMulti,
                  pairedSamples = pairedSamples
    )
    
    ## return the created object
    return(object)
}


## ---------
## Accessors
## ---------

#' @describeIn MirnaExperiment-class Access the results of miRNA differential
#' expression
#' @inheritParams deAccessors
#' @export
setMethod(
    "mirnaDE",
    "MirnaExperiment",
    function(object, onlySignificant, param, returnObject) {
        if (onlySignificant == TRUE &
            returnObject == FALSE &
            param == FALSE) {
            object@mirnaDE$data[object@mirnaDE$data$ID %in%
                                    object@mirnaDE$significant, ]
        } else if (onlySignificant == FALSE &
                   returnObject == FALSE &
                   param == FALSE) {
            object@mirnaDE$data
        } else if (param == TRUE &
                   returnObject == FALSE) {
            object@mirnaDE
        } else {
            object@mirnaDE$deObject
        }
    }
)

#' @describeIn MirnaExperiment-class Access the results of gene differential
#' expression
#' @inheritParams deAccessors
#' @export
setMethod(
    "geneDE",
    "MirnaExperiment",
    function(object, onlySignificant, param, returnObject) {
        if (onlySignificant == TRUE &
            returnObject == FALSE &
            param == FALSE) {
            object@geneDE$data[object@geneDE$data$ID %in%
                                   object@geneDE$significant, ]
        } else if (onlySignificant == FALSE &
                   returnObject == FALSE &
                   param == FALSE) {
            object@geneDE$data
        } else if (param == TRUE &
                   returnObject == FALSE) {
            object@geneDE
        } else {
            object@geneDE$deObject
        }
    }
)

#' @describeIn MirnaExperiment-class Access the names of differentially
#' expressed miRNAs
#' @inheritParams significantAccessors
#' @export
setMethod("significantMirnas", "MirnaExperiment", function(object) {
    object@mirnaDE$significant
})

#' @describeIn MirnaExperiment-class Access the names of differentially
#' expressed genes
#' @inheritParams significantAccessors
#' @export
setMethod("significantGenes", "MirnaExperiment", function(object) {
    object@geneDE$significant
})

#' @describeIn MirnaExperiment-class Check if the object derives from
#' sample-matched data
#' @inheritParams pairedSamples
#' @export
setMethod("pairedSamples", "MirnaExperiment", function(object) {
    object@pairedSamples
})

#' @describeIn MirnaExperiment-class Extract the miRNA-targets interactions
#' retrieved for the differentially expressed miRNAs
#' @inheritParams mirnaTargets
#' @export
setMethod("mirnaTargets", "MirnaExperiment", function(object) {
    object@targets
})

#' @describeIn MirnaExperiment-class Access the results of the integrative
#' miRNA-mRNA analysis
#' @inheritParams integration
#' @export
setMethod(
    "integration", "MirnaExperiment",
    function(object, param) {
        if (param == TRUE) {
            object@integration
        } else {
            object@integration$data
        }
    }
)


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

setReplaceMethod(
    "integration",
    "MirnaExperiment",
    function(object, value) {
        object@integration <- value
        validObject(object)
        object
    }
)



## ==========================================================================
## FunctionalEnrichment Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The `FunctionalEnrichment` class
#'
#' This class introduces the possibility to store the results of functional
#' enrichment analyses such as over-representation analysis (ORA), gene set
#' enrichment analysis (GSEA), and competitive gene set test accounting for
#' inter-gene correlation (CAMERA). The different slots contained in this class
#' are used to store enrichment results generated through the [enrichGenes()]
#' function.
#'
#' @slot data A `data.frame` object holding the output of enrichment analysis
#' @slot method The method used to perform functional enrichment analysis
#' (e.g. `Gene Set Enrichment Analysis (GSEA)`)
#' @slot organism The name of the organism under consideration (e.g.
#' `Homo sapiens`)
#' @slot database The name of the database used for the enrichment analysis
#' (e.g. `KEGG`)
#' @slot pCutoff A `numeric` value defining the threshold used for
#' statistical significance in the enrichment analysis (e.g. `0.05`)
#' @slot pAdjustment A `character` indicating the method used to correct
#' p-values for multiple testing (e.g. `fdr`)
#' @slot features A `character` vector containing the list of features used for
#' the enrichment
#' @slot statistic A `numeric` vector containing the statistic used to run
#' GSEA. This parameter is empty for ORA and CAMERA
#' @slot universe The background universe of features. Typically, this is equal
#' to the complete list of features assayed. This slot is NULL for GSEA
#' @slot geneSet The gene set used for the functional enrichment analysis. It
#' is a `list` object where each element contains the list of genes belonging
#' to a specific pathway.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name FunctionalEnrichment-class
#' @aliases FunctionalEnrichment
#' @docType class
#' @export
#' @import methods
setClass(
    "FunctionalEnrichment",
    representation(
        data = "data.frame",
        method = "character",
        organism = "character",
        database = "character",
        pCutoff = "numeric",
        pAdjustment = "character",
        features = "character",
        statistic = "numeric",
        universe = "character",
        geneSet = "list"
    )
)


## --------
## Validity
## --------

setValidity("FunctionalEnrichment", function(object) {
    if (!is.data.frame(object@data)) {
        return(paste(
            "'data' slot must be a data.frame that stores the results of",
            "functional enrichment analyses. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@method)) {
        return(paste(
            "'method' slot must be a character object that specifies the",
            "functonal enrichment method used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@organism)) {
        return(paste(
            "'organism' slot must be a character object that specifies",
            "the organism used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@database)) {
        return(paste(
            "'database' slot must be a character object that specifies",
            "the database used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.numeric(object@pCutoff)) {
        return(paste(
            "'pCutoff' slot must be a numeric object that specifies",
            "the p-value cutoff used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@pAdjustment) |
               !object@pAdjustment %in% p.adjust.methods) {
        return(paste(
            "'pAdjustment' slot must be a character object that specifies",
            "the p-value correction method used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@features)) {
        return(paste(
            "'features' slot must be a character object containing the",
            "features used for functional enrichment analysis. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.numeric(object@statistic)) {
        return(paste(
            "'statistic' slot must be a numeric object that specifies",
            "the metric used for GSEA. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else if (!is.character(object@universe)) {
        return(paste(
            "'universe' slot must be a character object that contains",
            "the complete list of background genes used. Please see",
            "?FunctionalEnrichment-class"
        ))
    } else {
        return(TRUE)
    }
})


## ---------
## Accessors
## ---------

#' @describeIn FunctionalEnrichment-class Access the `data` slot to take a
#' closer look at all the enriched terms of an enrichment analysis
#' @inheritParams enrichmentResults
#' @export
setMethod("enrichmentResults", "FunctionalEnrichment", function(object) {
    object@data
})

#' @describeIn FunctionalEnrichment-class See the database used for the
#' functional enrichment
#' @inheritParams enrichmentDatabase
#' @export
setMethod("enrichmentDatabase", "FunctionalEnrichment", function(object) {
    object@database
})

#' @describeIn FunctionalEnrichment-class Visualize the approach used for the
#' functional enrichment analysis
#' @inheritParams enrichmentMethod
#' @export
setMethod("enrichmentMethod", "FunctionalEnrichment", function(object) {
    object@method
})

#' @describeIn FunctionalEnrichment-class Access the `geneSet` slot to see
#' the collection of gene sets used for GSEA
#' @inheritParams geneSet
#' @export
setMethod("geneSet", "FunctionalEnrichment", function(object) {
    object@geneSet
})

#' @describeIn FunctionalEnrichment-class View the ranking metric used for GSEA
#' @inheritParams enrichmentMetric
#' @export
setMethod("enrichmentMetric", "FunctionalEnrichment", function(object) {
    object@statistic
})

#' @describeIn FunctionalEnrichment-class View the names of the pre-ranked
#' features used for GSEA
#' @inheritParams enrichedFeatures
#' @export
setMethod("enrichedFeatures", "FunctionalEnrichment", function(object) {
    object@features
})


## -------
## Setters
## -------

setReplaceMethod(
    "enrichmentResults", "FunctionalEnrichment",
    function(object, value) {
        object@data <- value
        validObject(object)
        object
    }
)

setReplaceMethod(
    "enrichmentDatabase", "FunctionalEnrichment",
    function(object, value) {
        object@database <- value
        validObject(object)
        object
    }
)



## ==========================================================================
## IntegrativePathwayAnalysis Class
## ==========================================================================


## ----------------
## Class definition
## ----------------


#' The `IntegrativePathwayAnalysis` class
#'
#' This class stores the output of integrative multi-omic pathway analyses.
#' In particular, the slots of this class are suitable to represent the results
#' of topologically-aware integrative pathway analysis (TAIPA) returned from
#' the [topologicalAnalysis()] function.
#'
#' @details
#'
#' ## Analysis results
#'
#' The `data` slot of this class consists in a `data.frame` object with six
#' columns, namely:
#'
#' * `pathway`, which indicates the name of the biological network;
#' * `coverage`, which specifies the fraction of nodes with expression
#' measurement available;
#' * `score`, which expresses the score of each individual pathway;
#' * `normalized.score`, which indicates the pathway scores after standardizing
#' the values for the null distribution computed through permutations;
#' * `P.Val`, the resulting p-value of each pathway;
#' * `adj.P.Val`, the p-value adjusted for multiple testing.
#'
#' ## Organisms and databases
#'
#' The `organism` and `database` slots specify the organism in study and the
#' database used for retrieving biological interactions, respectively. In
#' particular, the [topologicalAnalysis()] function supports `KEGG`,
#' `WikiPathways`, and `Reactome` databases. Regarding organisms, the
#' [supportedOrganisms()] function can be used to retrieve the available
#' species for each database.
#'
#' ## Statistical significance of the permutation test
#'
#' `pCutoff` and `pAdjustment` slots refer to the cutoff used for the analysis.
#' `pCutoff` is the threshold used for defining statistically significant
#' pathways, whereas `pAdjustment` refers to the multiple testing correction
#' method used.
#'
#' Furthermore, since the statistical significance of each pathway is defined
#' on the basis of a permutation test, the number of permutations is also
#' specified in the `nPerm` slot.
#'
#' ## Augmented pathways
#'
#' The `pathways` slot contains a `list` with weighted `graph` objects, each
#' representing a biological pathway. These networks have been enlarged by
#' adding the observed miRNA-mRNA interactions. Each network has been
#' processed so that the weight of each edge is +1 for activation interactions,
#' and -1 for repression interactions, such as those occurring between miRNAs
#' and mRNAs.
#'
#' ## Differential expression results for both miRNAs and genes
#'
#' The expression variation of all miRNAs and genes measured in the study is
#' stored in the `expression` slot. In particular, this slot consists of a
#' `data.frame` object with different information, including log2 fold changes,
#' node weights and p-values.
#'
#' ## Minimum percentage of measured features
#'
#' The `minPc` slot indicates the minimum percentage of miRNAs/mRNAs above
#' which pathways have been considered for the integrative analysis. This is
#' needed because often, when differential expression analysis is performed,
#' lowly expressed features are removed. Therefore, some pathways might result
#' significantly affected even if only 1% of nodes is perturbed.
#'
#' @slot data A `data.frame` object that contains the results of the
#' integrative pathway analysis. See the *details* section for further details
#' @slot method The method used for the analysis
#' @slot organism The name of the organism under consideration (e.g.
#' `Homo sapiens`)
#' @slot database The name of the database used for retrieving biological
#' pathways (e.g. `KEGG`)
#' @slot pCutoff A `numeric` value defining the threshold used for
#' statistical significance (e.g. `0.05`)
#' @slot pAdjustment A `character` indicating the method used to correct
#' p-values for multiple testing (e.g. `fdr`)
#' @slot pathways A `list` of `graph` objects containing the biological
#' networks retrieved from `database`, and augmented with
#' miRNA-mRNA interactions
#' @slot expression A `data.frame` object containing differential expression
#' results for both miRNAs and genes
#' @slot minPc The minimum percentage of measured features that a pathway must
#' have for being considered in the analysis
#' @slot nPerm The number of permutation used for assessing the statistical
#' significance of each pathway
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name IntegrativePathwayAnalysis-class
#' @aliases IntegrativePathwayAnalysis
#' @docType class
#' @export
#' @import methods
setClass(
    "IntegrativePathwayAnalysis",
    representation(
        data = "data.frame",
        method = "character",
        organism = "character",
        database = "character",
        pCutoff = "numeric",
        pAdjustment = "character",
        pathways = "list",
        expression = "data.frame",
        minPc = "numeric",
        nPerm = "numeric"
    )
)


## --------
## Validity
## --------

setValidity("IntegrativePathwayAnalysis", function(object) {
    if (!is.data.frame(object@data)) {
        return(paste(
            "'data' slot must be a data.frame that stores the results of",
            "an integrative pathway analysis. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.character(object@method)) {
        return(paste(
            "'method' slot must be a character object that describes the",
            "approach used for the integrative pathway analysis. Please",
            "see ?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.character(object@organism)) {
        return(paste(
            "'organism' slot must be a character object that specifies",
            "the organism used. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.character(object@database)) {
        return(paste(
            "'database' slot must be a character object that specifies",
            "the database used. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.numeric(object@pCutoff)) {
        return(paste(
            "'pCutoff' slot must be a numeric object that specifies",
            "the p-value cutoff used. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.character(object@pAdjustment) |
               !object@pAdjustment %in% c(p.adjust.methods, "max-T")) {
        return(paste(
            "'pAdjustment' slot must be a character object that specifies",
            "the p-value correction method used. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.list(object@pathways)) {
        return(paste(
            "'pathways' slot must be a list object containing graph",
            "objects with augmented biological networks. Please",
            "see ?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.data.frame(object@expression)) {
        return(paste(
            "'expression' slot must be a data.frame object containing",
            "differential expression results of miRNAs and genes. Please",
            "see ?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.numeric(object@minPc)) {
        return(paste(
            "'minPc' slot must be a numeric object that specifies",
            "the minimum percentage of measured nodes required. Please",
            "see ?IntegrativePathwayAnalysis-class"
        ))
    } else if (!is.numeric(object@nPerm)) {
        return(paste(
            "'nPerm' slot must be a numeric object that defines",
            "the number of permutations used. Please see",
            "?IntegrativePathwayAnalysis-class"
        ))
    } else {
        return(TRUE)
    }
})


## ---------
## Accessors
## ---------

#' @describeIn IntegrativePathwayAnalysis-class Access the results of
#' integrative miRNA-mRNA pathway analysis
#' @inheritParams integratedPathways
#' @export
setMethod(
    "integratedPathways",
    "IntegrativePathwayAnalysis",
    function(object) {
        object@data
    }
)

#' @describeIn IntegrativePathwayAnalysis-class View the database used for
#' the integrative pathway analysis
#' @inheritParams integrationDatabase
#' @export
setMethod(
    "integrationDatabase",
    "IntegrativePathwayAnalysis",
    function(object) {
        object@database
    }
)

#' @describeIn IntegrativePathwayAnalysis-class Extract the list of biological
#' networks augmented with miRNA-mRNA interactions
#' @inheritParams augmentedPathways
#' @export
setMethod(
    "augmentedPathways",
    "IntegrativePathwayAnalysis",
    function(object) {
        object@pathways
    }
)


## -------
## Setters
## -------

setReplaceMethod(
    "integratedPathways",
    "IntegrativePathwayAnalysis",
    function(object, value) {
        object@data <- value
        validObject(object)
        object
    }
)
