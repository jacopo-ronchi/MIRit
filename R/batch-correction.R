#' Correct for batch effects in miRNA and gene expression measurements
#'
#' This function allows to remove unwanted batch effects from miRNA and gene
#' expression matrices. In particular, this function fits a linear model to
#' miRNA/gene expression levels, and then removes the variability caused by
#' batch effects. Furthermore, a weighted surrogate variable analysis (WSVA)
#' can also be included to remove the effects due to surrogate variables.
#' If batch effects are present, it is crucial to remove them with this
#' function before moving to correlation analysis.
#'
#' @details
#' Batch effects consist in unwanted sources of technical variation that
#' confound expression variability and limit downstream analyses. Since the
#' reliability of biological conclusions of integrative miRNA-mRNA analyses
#' depends on the association between miRNA and gene expression levels, it is
#' pivotal to ensure that expression measurements are not affected by technical
#' variations. In this regard, if batch effects are noticed in the data, the
#' user should run this function before using the [mirnaIntegration()]
#' function to perform a correlation analysis.
#'
#' Usually, given a [`MirnaExperiment`][MirnaExperiment-class] object, the user
#' should specify:
#'
#' * the `assay` from which we want to remove batch effects (one between
#' `genes` and `microRNA`);
#' * the `batch` variable, which is a variable that defines the different
#' batches;
#' * the `batch2` variable, which can be included to correct for a second
#' series of batches that have additive effects to those specified in `batch`;
#' * the `covariates` variables, which allows correction for one or more
#' continuous numeric effects.
#'
#' In particular, `batch` and `batch2` could be provided as the names of
#' covariates included in the `colData` of a
#' [`MirnaExperiment`][MirnaExperiment-class] object. Alternatively, they can
#' be `character`/`factor` objects that declare batch memberships.
#' Similarly, `covariates` can be supplied as a vector containing the names
#' of numeric variables listed in the `colData` of
#' [`MirnaExperiment`][MirnaExperiment-class] objects, or they can be provided
#' as a simple `matrix`.
#'
#' Additionally, the influence of unknown sources of technical variation can
#' be removed by including surrogate variables estimated through WSVA. To do
#' so, we can set `includeWsva` to TRUE, and then we can specify the number of
#' surrogate variables to use through the `n.sv` parameter. Further, the
#' surrogate variables can be tuned to the more variable genes by setting
#' `weight.by.sd` to TRUE.
#'
#' Please note that we only recommend to remove batch effects directly from
#' expression measurements prior to correlation analysis. This function
#' can't be used to remove batch effects before differential expression
#' analysis, because for that purpose, it is better to include batch variables
#' in the linear model. In this way, we do not underestimate the residual
#' degrees of freedom, so that the calculated standard errors, t-statistics
#' and p-values are not overoptimistic.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param assay The expression matrix to correct. It must be one of `genes`
#' and `microRNA`
#' @param batch It must be the name of a variable present in the `colData` of
#' a [`MirnaExperiment`][MirnaExperiment-class] object (eg. "disease"), or,
#' alternatively, it must be a `character`/`factor` object that defines batch
#' memberships. See the **details** section for additional information
#' @param batch2 It must be the name of a variable present in the `colData` of
#' a [`MirnaExperiment`][MirnaExperiment-class] object (eg. "disease"), or,
#' alternatively, it must be a `character`/`factor` object that defines
#' another series of batches that have additive effects to those specified
#' in `batch`. See the **details** section for additional information
#' @param covariates Additional numeric covariates that we want to correct for.
#' It must be a `character` vector containing the names of numeric variables
#' present in the `colData` of a [`MirnaExperiment`][MirnaExperiment-class]
#' object (eg. `c("age", "RIN", "quantity")`), or, alternatively, it must be a
#' simple `matrix` object. See the **details** section for additional
#' information
#' @param includeWsva Logical, whether to correct for surrogate variables or
#' not. Default is FALSE
#' @param n.sv The number of surrogate variables to estimate
#' @param weight.by.sd Logical, whether to specifically tune the surrogate
#' variables to the more variable genes or not. Default is TRUE
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing batch
#' effect-corrected expression matrices.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # correct batch effects due to the patient from miRNA expression matrix
#' obj <- batchCorrection(obj, "microRNA", batch = "patient")
#'
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies.” Nucleic Acids Research, 43(7), e47. \url{doi:10.1093/nar/gkv007}.
#'
#' @note
#' To estimate surrogate variables and to remove batch effects from expression
#' data, MIRit uses the [limma::wsva()] and [limma::removeBatchEffect()]
#' functions, respectively.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
batchCorrection <- function(mirnaObj,
                            assay,
                            batch = NULL,
                            batch2 = NULL,
                            covariates = NULL,
                            includeWsva = FALSE,
                            n.sv = 1L,
                            weight.by.sd = TRUE) {
    ## check inputs
    if (!is(mirnaObj, "MirnaExperiment")) {
        stop("'mirnaObj' should be of class MirnaExperiment! ",
             "See ?MirnaExperiment",
             call. = FALSE
        )
    }
    if (length(assay) != 1 |
        !assay %in% c("microRNA", "genes")) {
        stop("'assay' must be either 'microRNA' or 'genes'!", call. = FALSE)
    }
    if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0 &
        assay == "microRNA") {
        stop("MiRNA differential expression results are not present in ",
             "'mirnaObj'. Please, use 'performMirnaDE()' before using ",
             "this function. See ?performMirnaDE",
             call. = FALSE
        )
    }
    if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0 &
        assay == "genes") {
        stop("Gene differential expression results are not present in ",
             "'mirnaObj'. Please, use 'performGeneDE()' before using ",
             "this function. See ?performGeneDE",
             call. = FALSE
        )
    }
    if (!is.null(batch)) {
        if (length(batch) == 1) {
            if (!is.character(batch) |
                !(batch %in% colnames(colData(mirnaObj)) &
                  !batch %in% c("primary", "mirnaCol", "geneCol"))) {
                stop("'batch' must be a batch variable present in the ",
                     "metadata (colData) of a MirnaExperiment object; or, ",
                     "alternatively, it must be a character/factor object ",
                     "that specifies the different batches.",
                     call. = FALSE
                )
            }
        } else {
            if ((!is.character(batch) & !is.factor(batch)) |
                length(batch) != ncol(mirnaObj[[assay]])) {
                stop("'batch' must be a batch variable present in the ",
                     "metadata (colData) of a MirnaExperiment object; or, ",
                     "alternatively, it must be a character/factor object ",
                     "that specifies the different batches.",
                     call. = FALSE
                )
            }
        }
    }
    if (!is.null(batch2)) {
        if (length(batch2) == 1) {
            if (!is.character(batch2) |
                !(batch2 %in% colnames(colData(mirnaObj)) &
                  !batch2 %in% c("primary", "mirnaCol", "geneCol"))) {
                stop("'batch2' must be a batch variable present in the ",
                     "metadata (colData) of a MirnaExperiment object; or, ",
                     "alternatively, it must be a character/factor object ",
                     "that specifies the different batches.",
                     call. = FALSE
                )
            }
        } else {
            if ((!is.character(batch2) & !is.factor(batch2)) |
                length(batch2) != ncol(mirnaObj[[assay]])) {
                stop("'batch2' must be a batch variable present in the ",
                     "metadata (colData) of a MirnaExperiment object; or, ",
                     "alternatively, it must be a character/factor object ",
                     "that specifies the different batches.",
                     call. = FALSE
                )
            }
        }
    }
    if (!is.null(covariates)) {
        if (is.character(covariates)) {
            if (any(!covariates %in% colnames(colData(mirnaObj)))) {
                stop("'covariates' must be either a vector of column names ",
                     "present in the metadata (colData) of a MirnaExperiment ",
                     "object; or, alternatively, it must be a numeric vector ",
                     "or matrix object that specifies the covariates to ",
                     "correct for.",
                     call. = FALSE
                )
            }
        } else if (is.numeric(covariates)) {
            if (is.matrix(covariates)) {
                if (nrow(covariates) != ncol(mirnaObj[[assay]])) {
                    stop("'covariates' must be either a vector of column ",
                         "names present in the metadata (colData) of a ",
                         "MirnaExperiment object; or, alternatively, it must ",
                         "be a numeric vector or matrix object that ",
                         "specifies the covariates to correct for.",
                         call. = FALSE
                    )
                }
            } else {
                if (length(covariates) != ncol(mirnaObj[[assay]])) {
                    stop("'covariates' must be either a vector of column ",
                         "names present in the metadata (colData) of a ",
                         "MirnaExperiment object; or, alternatively, it must ",
                         "be a numeric vector or matrix object that ",
                         "specifies the covariates to correct for.",
                         call. = FALSE
                    )
                }
            }
        } else {
            stop("'covariates' must be either a vector of column names ",
                 "present in the metadata (colData) of a MirnaExperiment ",
                 "object; or, alternatively, it must be a numeric vector ",
                 "or matrix object that specifies the covariates to ",
                 "correct for.",
                 call. = FALSE
            )
        }
    }
    if (!is.logical(includeWsva) |
        length(includeWsva) != 1) {
        stop("'includeWsva' must be logical (TRUE/FALSE)!", call. = FALSE)
    }
    if (!is.numeric(n.sv) |
        length(n.sv) != 1 |
        n.sv < 0 |
        !n.sv %% 1 == 0) {
        stop("'n.sv' must be a non-neagtive number! (default is 1)",
             call. = FALSE
        )
    }
    if (!is.logical(weight.by.sd) |
        length(weight.by.sd) != 1) {
        stop("'weight.by.sd' must be logical (TRUE/FALSE)!", call. = FALSE)
    }
    
    ## verify that batch correction has not been performed yet
    if (!is.null(mirnaObj@metadata$uncorrectedMatrices[[assay]])) {
        warning("Batch effect correction has already been performed for ",
                "this assay! The uncorrected matrix has been restored ",
                "before overwriting the batch-corrected matrix.",
                call. = FALSE
        )
        mirnaObj[[assay]] <- mirnaObj@metadata$uncorrectedMatrices[[assay]]
    } else {
        ## move normalized counts to metadata
        mirnaObj@metadata$uncorrectedMatrices[[assay]] <- mirnaObj[[assay]]
    }
    
    ## extract expression matrix
    data <- mirnaObj[[assay]]
    
    ## extract the parameters used for differential expression analysis
    if (assay == "microRNA") {
        deParam <- mirnaDE(mirnaObj, param = TRUE)
    } else if (assay == "genes") {
        deParam <- geneDE(mirnaObj, param = TRUE)
    }
    design <- deParam$design
    group <- deParam$group
    
    ## extract the covariates to correct for
    sm <- sampleMap(mirnaObj)
    sm <- sm[sm$assay == assay, ]
    rownames(sm) <- sm$colname
    sm <- sm[colnames(mirnaObj[[assay]]), ]
    if (is.character(covariates)) {
        smCov <- sm[, covariates]
        allNumeric <- vapply(colnames(smCov),
                             function(x) is.numeric(smCov[, x]),
                             FUN.VALUE = logical(1)
        )
        if (any(allNumeric == FALSE)) {
            stop("'covariates' must only contain numeric variables! ",
                 "See ?batchCorrection",
                 call. = FALSE
            )
        }
        covariates <- as.matrix(smCov)
    }
    
    ## perform surrogate variable analysis
    if (includeWsva == TRUE) {
        des <- model.matrix(design, data = sm)
        ws <- limma::wsva(data,
                          design = des,
                          n.sv = n.sv, weight.by.sd = weight.by.sd
        )
        if (!is.null(covariates)) {
            covariates <- cbind(covariates, ws)
        } else {
            covariates <- ws
        }
    }
    
    ## define the batch and group effects
    if (!is.null(batch)) {
        batch <- sm[, batch]
    }
    if (!is.null(batch2)) {
        batch2 <- sm[, batch2]
    }
    
    ## remove batch effects from expression matrix
    data <- limma::removeBatchEffect(data,
                                     batch = batch,
                                     batch2 = batch2,
                                     covariates = covariates,
                                     group = sm[, group]
    )
    
    ## set expression levels to the batch-corrected ones
    mirnaObj[[assay]] <- data
    
    ## return the MirnaExperiment object
    return(mirnaObj)
}
