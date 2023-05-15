#' Perform differential expression analysis
#' 
#' `performMirnaDE()` and `performGeneDE()` are two functions provided by MIRit
#' to conduct miRNA and gene differential expression analysis, respectively.
#' In particular, these functions allow the user to compute differential
#' expression through different methods, namely `edgeR`, `DESeq2`, `limma-voom`
#' and `limma`. Data deriving from NGS experiments and microarray technology
#' are all suitable for these functions. For precise indications about how to
#' use these functions, please refer to the *details* section.
#' 
#' @details
#' When performing differential expression for NGS experiments, count matrices
#' are detected and `method` parameter must be one of `edgeR`, `DESeq2`,
#' and `voom`. On the other hand, when dealing with microarray studies, only
#' `limma` can be used.
#' 
#' To calculate differential expression, MIRit must be informed about the
#' variable of interest and the desired contrast. In particular, the `group`
#' parameter must be the name of a variable present in the metadata (colData)
#' of a [`MirnaExperiment`][MirnaExperiment-class] object, which specifies the
#' variable used to compute differential expression analysis, between the groups
#' indicated in `contrast`. Specifically, `contrast` must be a character vector
#' that defines the levels to compare separated by a dash. For example, if we
#' have a variable named 'condition', with two levels, namely 'disease' and
#' 'healthy', we can identify differentially expressed genes in 'disease'
#' samples compared to 'healthy' subjects by specifying: `group = 'condition'`
#' and `contrast = 'disease-healthy'`. Furthermore, the user needs to specify
#' the model to fit expression values. To do so, the user has to state the
#' model formula in the `design` parameter. Please note that for a correct
#' inner working of these functions, the `group` variable of interest must be
#' the *first* variable in model formula. Moreover, the user can include in the
#' design any other sources of variation by specifying covariates that will be
#' taken into account. For instance, if we want to compare 'disease' subjects
#' against 'healthy' individuals, without the influence of sex differences,
#' we may specify `design = ~ condition + sex`, where 'sex' is also a
#' variable present in the metadata (colData) of `mirnaObj`.
#' 
#' Notably, for all the methods available, the user can supply additional
#' arguments to the functions implemented in `edgeR`, `DESeq2` and `limma`.
#' Therefore, the user has finer control over how the differential expression
#' analysis is performed. In this regard, for microarray studies, the user
#' may opt to include weighted surrogate variable analysis (WSVA) to correct
#' for unknown sources of variation (`useWsva = TRUE`). Moreover, for
#' microarray data, the `arrayWeights()` function in `limma` can be used to
#' assess differential expression with respect to array qualities. Additionally,
#' when using `limma-voom`, the user may estimate voom trasformation with or
#' without quality weights (by specifying `useVoomWithQualityWeights = TRUE`).
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param group The variable of interest for differential expression analysis.
#' It must be the column name of a variable present in the metadata (colData)
#' of a [`MirnaExperiment`][MirnaExperiment-class] object. See the *details*
#' section for additional information
#' @param contrast A `character` object that specifies the groups to be
#' compared during differential expression analysis, separated by a dash
#' (e.g. 'disease-healthy'). Note that reference group must be the last one,
#' for additional information see the *details* section
#' @param design An R `formula` that indicates the model to fit. It must
#' include the variable of interest (`group`) together with eventual
#' covariates (e.g. '~ 0 + disease + sex'). Please note that `group` variable
#' must be the first one. See the *details* section for additional information
#' @param method The statistical package used to compute differential
#' expression. For NGS experiments, it must be one of `edgeR` (default),
#' `DESeq2`, and `voom` (for limma-voom). Instead, for microarray data, only
#' `limma` can be used
#' @param logFC The minimum log2 fold change required to consider a gene as
#' differentially expressed. Default is 1, to retain only two-fold differences
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param filterByExpr.args A `list` object containing additional arguments
#' passed to [edgeR::filterByExpr()] function. It is used when `method` is set
#' to `edgeR` or `voom`
#' @param calcNormFactors.args A `list` object containing additional arguments
#' passed to [edgeR::calcNormFactors()] function. It is used when `method` is
#' set to `edgeR` or `voom`
#' @param estimateDisp.args A `list` object containing additional arguments
#' passed to [edgeR::estimateDisp()] function. It is used when `method` is
#' set to `edgeR`. Default is `list(robust = TRUE)` to use the robust parameter
#' @param glmQLFit.args A `list` object containing additional arguments
#' passed to [edgeR::glmQLFit()] function. It is used when `method` is
#' set to `edgeR`
#' @param glmQLFTest.args A `list` object containing additional arguments
#' passed to [edgeR::glmQLFTest()] function. It is used when `method` is
#' set to `edgeR`
#' @param DESeq.args A `list` object containing additional arguments
#' passed to [DESeq2::DESeq()] function. It is used when `method` is
#' set to `DESeq`
#' @param useVoomWithQualityWeights Logical, whether to use the
#' [limma::voomWithQualityWeights()] function or just the [limma::voom()]
#' function. It is used when `method` is set to `voom`. Default is TRUE
#' @param voom.args A `list` object containing additional arguments
#' passed to [limma::voom()] function or [limma::voomWithQualityWeights()]
#' function. It is used when `method` is set to `voom`
#' @param lmFit.args A `list` object containing additional arguments
#' passed to [limma::lmFit()] function. It is used when `method` is set
#' to `voom` or `limma`
#' @param eBayes.args A `list` object containing additional arguments
#' passed to [limma::eBayes()] function. It is used when `method` is set
#' to `voom` or `limma`
#' @param useArrayWeights Logical, whether to use the [limma::arrayWeights()]
#' function or not. It is used when `method` is set to `limma`. Default is TRUE
#' @param useWsva Logical, whether to use the [limma::wsva()] function or not.
#' It is used when `method` is set to `limma`. Default is FALSE
#' @param arrayWeights.args A `list` object containing additional arguments
#' passed to [limma::arrayWeights()] function. It is used when `method` is set
#' to `limma`
#' @param wsva.args A `list` object containing additional arguments
#' passed to [limma::wsva()] function. It is used when `method` is set
#' to `limma`
#' 
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing differential
#' expression results. To access these results, the user may run the
#' [mirnaDE()] and [geneDE()] functions for miRNAs and genes, respectively.
#' 
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform miRNA DE with edgeR
#' obj <- performMirnaDE(obj, group = "disease", contrast = "PTC-NTH",
#' design = ~ 0 + disease + patient, method = "edgeR")
#' 
#' # perform miRNA DE with DESeq2
#' obj <- performMirnaDE(obj, group = "disease", contrast = "PTC-NTH",
#' design = ~ 0 + disease + patient, method = "DESeq2")
#' 
#' # perform gene DE with limma-voom
#' obj <- performGeneDE(obj, group = "disease", contrast = "PTC-NTH",
#' design = ~ 0 + disease + patient, method = "voom")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @name deAnalysis
NULL





#' @describeIn deAnalysis Perform differential expression analysis for miRNAs
#' @export
performMirnaDE <- function(mirnaObj,
                           group,
                           contrast,
                           design,
                           method = "edgeR",
                           logFC = 1,
                           pCutoff = 0.05,
                           pAdjustment = "fdr",
                           filterByExpr.args = list(),
                           calcNormFactors.args = list(),
                           estimateDisp.args = list(robust = TRUE),
                           glmQLFit.args = list(),
                           glmQLFTest.args = list(),
                           DESeq.args = list(),
                           useVoomWithQualityWeights = TRUE,
                           voom.args = list(),
                           lmFit.args = list(),
                           eBayes.args = list(),
                           useArrayWeights = TRUE,
                           useWsva = FALSE,
                           wsva.args = list(),
                           arrayWeights.args = list()) {
  
  ## perform differential expression analyis for miRNAs
  mirnaObj <- performDE(assay = "miRNAs",
                        mirnaObj = mirnaObj,
                        group = group,
                        contrast = contrast,
                        design = design,
                        method = method,
                        logFC = logFC,
                        pCutoff = pCutoff,
                        pAdjustment = pAdjustment,
                        filterByExpr.args = filterByExpr.args,
                        calcNormFactors.args = calcNormFactors.args,
                        estimateDisp.args = estimateDisp.args,
                        glmQLFit.args = glmQLFit.args,
                        glmQLFTest.args = glmQLFTest.args,
                        DESeq.args = DESeq.args,
                        useVoomWithQualityWeights = useVoomWithQualityWeights,
                        voom.args = voom.args,
                        lmFit.args = lmFit.args,
                        eBayes.args = eBayes.args,
                        useArrayWeights = useArrayWeights,
                        useWsva = useWsva,
                        wsva.args = wsva.args,
                        arrayWeights.args = arrayWeights.args)
  
  ## return mirnaObj
  return(mirnaObj)
  
}





#' @describeIn deAnalysis Perform differential expression analysis for genes
#' @export
performGeneDE <- function(mirnaObj,
                          group,
                          contrast,
                          design,
                          method = "edgeR",
                          logFC = 1,
                          pCutoff = 0.05,
                          pAdjustment = "fdr",
                          filterByExpr.args = list(),
                          calcNormFactors.args = list(),
                          estimateDisp.args = list(robust = TRUE),
                          glmQLFit.args = list(),
                          glmQLFTest.args = list(),
                          DESeq.args = list(),
                          useVoomWithQualityWeights = TRUE,
                          voom.args = list(),
                          lmFit.args = list(),
                          eBayes.args = list(),
                          useArrayWeights = TRUE,
                          useWsva = FALSE,
                          wsva.args = list(),
                          arrayWeights.args = list()) {
  
  ## perform differential expression analyis for genes
  mirnaObj <- performDE(assay = "genes",
                        mirnaObj = mirnaObj,
                        group = group,
                        contrast = contrast,
                        design = design,
                        method = method,
                        logFC = logFC,
                        pCutoff = pCutoff,
                        pAdjustment = pAdjustment,
                        filterByExpr.args = filterByExpr.args,
                        calcNormFactors.args = calcNormFactors.args,
                        estimateDisp.args = estimateDisp.args,
                        glmQLFit.args = glmQLFit.args,
                        glmQLFTest.args = glmQLFTest.args,
                        DESeq.args = DESeq.args,
                        useVoomWithQualityWeights = useVoomWithQualityWeights,
                        voom.args = voom.args,
                        lmFit.args = lmFit.args,
                        eBayes.args = eBayes.args,
                        useArrayWeights = useArrayWeights,
                        useWsva = useWsva,
                        wsva.args = wsva.args,
                        arrayWeights.args = arrayWeights.args)
  
  ## return mirnaObj
  return(mirnaObj)
  
}





## internal function to perform differential expression analysis
performDE <- function(assay,
                      mirnaObj,
                      group,
                      contrast,
                      design,
                      method = "edgeR",
                      logFC = 1,
                      pCutoff = 0.05,
                      pAdjustment = "fdr",
                      filterByExpr.args = list(),
                      calcNormFactors.args = list(),
                      estimateDisp.args = list(robust = TRUE),
                      glmQLFit.args = list(),
                      glmQLFTest.args = list(),
                      DESeq.args = list(),
                      useVoomWithQualityWeights = TRUE,
                      voom.args = list(),
                      lmFit.args = list(),
                      eBayes.args = list(),
                      useArrayWeights = TRUE,
                      useWsva = FALSE,
                      wsva.args = list(),
                      arrayWeights.args = list()) {
  
  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!is.character(group) |
      length(group) != 1 |
      !(group %in% colnames(MultiAssayExperiment::colData(mirnaObj)) &
        !group %in% c("primary", "mirnaCol", "geneCol"))) {
    stop(paste("'group' must be the column name of a variable specified",
               "in the metadata (colData) of a MirnaExperiment object."),
         call. = FALSE)
  }
  if (!is.character(contrast) |
      length(contrast) != 1 |
      length(strsplit(contrast, "-")[[1]]) != 2 |
      !all(strsplit(contrast, "-")[[1]] %in%
           MultiAssayExperiment::colData(mirnaObj)[, group])) {
    stop(paste("'contrast' must be a character that specifies the groups",
               "for which you want to calculate differential expression",
               "(e.g. 'PTC-NTH'). For details, see ?performMirnaDE or",
               "?performGeneDE"),
         call. = FALSE)
  }
  if (!rlang::is_formula(design) |
      !all(labels(stats::terms(design)) %in%
           colnames(MultiAssayExperiment::colData(mirnaObj)) &
           !labels(stats::terms(design)) %in%
           c("primary", "mirnaCol", "geneCol"))) {
    stop(paste("'design' must be an R formula that specifies the variables",
               "in colData that will be used to model expression. For",
               "details, see ?performMirnaDE or ?performGeneDE"),
         call. = FALSE)
  }
  if (!is.character(method) |
      length(method) != 1 |
      !method %in% c("limma", "edgeR", "DESeq2", "voom")) {
    stop(paste("'method' must be  one of: 'limma', 'edgeR' (default),",
               "'DESeq2', 'voom'. For additional details, see ?performMirnaDE",
               "ord ?performGeneDE"),
         call. = FALSE)
  }
  if (!is.numeric(logFC) |
      length(logFC) != 1 |
      logFC < 0) {
    stop(paste("'logFC' must be a non-neagtive number that specifies the",
               "minimum absolute significant fold change (default is 1)"),
         call. = FALSE)
  }
  if (!is.numeric(pCutoff) |
      length(pCutoff) != 1 |
      pCutoff > 1 |
      pCutoff < 0) {
    stop("'pCutoff' must be a number between 0 and 1! (default is 0.05)",
         call. = FALSE)
  }
  if (!is.character(pAdjustment) |
      length(pAdjustment) != 1 |
      !pAdjustment %in% c("none", "fdr", "bonferroni", "BY", "hochberg",
                          "holm", "hommel", "BH")) {
    stop(paste("'pAdjustment' must be  one of: 'none', 'fdr' (default),",
               "'BH' (same as 'fdr'), 'bonferroni', 'BY', 'hochberg',",
               "'holm', 'hommel'"),
         call. = FALSE)
  }
  if (!is.logical(useVoomWithQualityWeights) |
      length(useVoomWithQualityWeights) != 1) {
    stop("'useVoomWithQualityWeights' must be logical (TRUE/FALSE)!",
         call. = FALSE)
  }
  if (!is.logical(useArrayWeights) |
      length(useArrayWeights) != 1) {
    stop("'useArrayWeights' must be logical (TRUE/FALSE)!",
         call. = FALSE)
  }
  if (!is.logical(useWsva) |
      length(useWsva) != 1) {
    stop("'useWsva' must be logical (TRUE/FALSE)!",
         call. = FALSE)
  }
  if (!is.list(filterByExpr.args) |
      !is.list(calcNormFactors.args) |
      !is.list(estimateDisp.args) |
      !is.list(glmQLFit.args) |
      !is.list(glmQLFTest.args) |
      !is.list(DESeq.args) |
      !is.list(voom.args) |
      !is.list(lmFit.args) |
      !is.list(eBayes.args) |
      !is.list(wsva.args) |
      !is.list(arrayWeights.args)) {
    stop(paste("Additional arguments passed to the limma, edgeR and DESeq2",
               "functions must be passed as lists! See ?performMirnaDE",
               "or ?performGeneDE"),
         call. = FALSE)
  }
  
  ## define assay name and any previous DE object
  if (assay == "miRNAs") {
    assayName <- "microRNA"
    assayFunc <- "'mirnaDE()'"
    featCol <- "mirnaCol"
    prevDe <- mirnaDE(mirnaObj, returnObject = TRUE)
  } else if (assay == "genes") {
    assayName <- "genes"
    assayFunc <- "'geneDE()'"
    featCol <- "geneCol"
    prevDe <- geneDE(mirnaObj, returnObject = TRUE)
  }
  
  ## check if differential expression has already been carried out
  if (!is.null(prevDe)) {
    
    ## set expression back to counts
    mirnaObj[[assayName]] <-
      MultiAssayExperiment::metadata(mirnaObj)[["oldCounts"]][[assayName]]
    
  } else {
    
    ## move raw count matrices to metadata slot
    oldCounts <- list(mirnaObj[[assayName]])
    names(oldCounts) <- assayName
    MultiAssayExperiment::metadata(mirnaObj) <- list(oldCounts = oldCounts)
    
  }
  
  ## extract feature expression
  featExpr <- mirnaObj[[assayName]]
  
  ## extract sample metadata
  samplesMetadata <- MultiAssayExperiment::colData(mirnaObj)
  meta <- samplesMetadata[which(!is.na(samplesMetadata[, featCol])), ]
  
  ## determine if data derive from RNA-Seq or microarray experiments
  if (!all(featExpr%%1 == 0) & method != "limma") {
    warning(paste(method, "is not suitable for microarray experiments and",
                  "requires count data. Instead, 'limma' will be used",
                  "to assess differential expression..."), call. = FALSE)
    method <- "limma"
  } else if(all(featExpr%%1 == 0) & method == "limma") {
    warning(paste(method, "is not suitable for NGS experiments!",
                  "Instead, 'limma-voom' will be used",
                  "to assess differential expression..."), call. = FALSE)
    method <- "voom"
  } else {
    message(paste("Performing differential expression analysis with ",
                  method, "...", sep = ""))
  }
  
  ## perform differential expression with the appropriate method
  if (method == "edgeR") {
    deList <- edgeR.DE(counts = featExpr,
                       group = group,
                       contrast = contrast,
                       meta = meta,
                       design = design,
                       logFC = logFC,
                       pCutoff = pCutoff,
                       pAdjustment = pAdjustment,
                       filterByExpr.args = filterByExpr.args,
                       calcNormFactors.args = calcNormFactors.args,
                       estimateDisp.args = estimateDisp.args,
                       glmQLFit.args = glmQLFit.args,
                       glmQLFTest.args = glmQLFTest.args)
  } else if (method == "DESeq2") {
    deList <- DESeq2.DE(counts = featExpr,
                        group = group,
                        contrast = contrast,
                        meta = meta,
                        design = design,
                        logFC = logFC,
                        pCutoff = pCutoff,
                        pAdjustment = pAdjustment,
                        DESeq.args = DESeq.args)
  } else if (method == "voom") {
    deList <- voom.DE(counts = featExpr,
                      group = group,
                      contrast = contrast,
                      meta = meta,
                      design = design,
                      logFC = logFC,
                      pCutoff = pCutoff,
                      pAdjustment = pAdjustment,
                      useVoomWithQualityWeights = useVoomWithQualityWeights,
                      filterByExpr.args = filterByExpr.args,
                      calcNormFactors.args = calcNormFactors.args,
                      voom.args = voom.args,
                      lmFit.args = lmFit.args,
                      eBayes.args = eBayes.args)
  } else if (method == "limma") {
    deList <- limma.DE(expr = featExpr,
                       group = group,
                       contrast = contrast,
                       meta = meta,
                       design = design,
                       logFC = logFC,
                       pCutoff = pCutoff,
                       pAdjustment = pAdjustment,
                       useArrayWeights = useArrayWeights,
                       useWsva = useWsva,
                       wsva.args = wsva.args,
                       arrayWeights.args = arrayWeights.args,
                       lmFit.args = lmFit.args,
                       eBayes.args = eBayes.args)
  }
  
  ## change count matrices in mirnaObj with normalized expression matrices
  mirnaObj[[assayName]] <- deList[["normExpr"]]
  
  ## add differential expression results to mirnaObj
  if (assay == "miRNAs") {
    mirnaDE(mirnaObj) <- deList[seq(4)]
  } else if (assay == "genes") {
    geneDE(mirnaObj) <- deList[seq(4)]
  }
  
  ## inform the user about differential expression results
  message(paste("Differential expression analysis reported ",
                length(deList$significant), " significant ", assay,
                " with p < ", pCutoff, " (correction: ", pAdjustment, ").",
                " You can use the ", assayFunc, " function to access results.",
                sep = ""))
  
  ## return the object after differential expression
  return(mirnaObj)
  
}





## perform differential expression with edgeR
edgeR.DE <- function(counts,
                     group,
                     contrast,
                     meta,
                     design,
                     logFC,
                     pCutoff,
                     pAdjustment,
                     filterByExpr.args,
                     calcNormFactors.args,
                     estimateDisp.args,
                     glmQLFit.args,
                     glmQLFTest.args) {
  
  ## identify numerator level and reference level
  contrast <- strsplit(contrast, "-")[[1]]
  
  ## convert group variable to factor with specified reference level
  meta[, group] <- stats::relevel(factor(meta[, group]), ref = contrast[2])
  
  ## create edgeR object from counts
  features <- edgeR::DGEList(counts = counts,
                             group = meta[, group],
                             samples = meta)
  
  ## filter features based on expression
  keep <- do.call(edgeR::filterByExpr, c(list(features),
                                         filterByExpr.args))
  features <- features[keep, , keep.lib.sizes = FALSE]
  
  ## normalize feature counts (default with TMM)
  features <- do.call(edgeR::calcNormFactors,
                      c(list(features), calcNormFactors.args))
  
  ## design the model
  des <- stats::model.matrix(design, data = meta)
  
  ## estimate dispersion and fit negative binomial distribution
  features <- do.call(edgeR::estimateDisp,
                      c(list(features, des), estimateDisp.args))
  fit <- do.call(edgeR::glmQLFit, c(list(features, des), glmQLFit.args))
  
  ## determine if the supplied model has intercept
  intercept <- attributes(terms(design))["intercept"] == 1
  
  ## identify the comparison for DE analysis
  if (intercept == FALSE) {
    
    ## build the contrast of interest
    contrast <- paste(group, contrast, sep = "")
    contrast <- paste(contrast, collapse = "-")
    
    ## fit the contrast
    con <- limma::makeContrasts(contrasts = contrast,
                                levels = des)
    
    ## set contrast parameter for DE
    comparison <- list(contrast = con)
    
  } else {
    
    ## set the coefficient name for the appropriate comparison
    cf <- contrast[1]
    cf <- paste(group, cf, sep = "")
    
    ## set contrast parameter for DE
    comparison <- list(coef = cf)
    
  }
  
  ## perform differential expression
  de <- do.call(edgeR::glmQLFTest,
                c(list(fit), comparison, glmQLFTest.args))
  
  ## extract normalized expression matrix
  normExpr <- edgeR::cpm(features, normalized.lib.sizes = TRUE, log = TRUE)
  
  ## extract differential expression results
  deRes <- as.data.frame(edgeR::topTags(de, n = Inf,
                                        adjust.method = pAdjustment))
  
  ## add 'ID' column to differentially expressed results
  deRes$ID <- rownames(deRes)
  
  ## format differential expression table as required by MIRit
  deRes <- identifyColNames(deRes)
  
  ## select significant features
  sig <- rownames(deRes[abs(deRes$logFC) > logFC &
                          deRes$adj.P.Val < pCutoff, ])
  
  ## define the parameters used
  met <- paste("edgeR (p < ", pCutoff, ", correction: ",
               pAdjustment, ")", sep = "")
  
  ## create a list with DE results
  deList <- list(data = deRes,
                 significant = sig,
                 method = met,
                 deObject = features,
                 normExpr = normExpr)
  
  ## return differential expression results
  return(deList)
  
}





## perform differential expression with DESeq2
DESeq2.DE <- function(counts,
                      group,
                      contrast,
                      meta,
                      design,
                      logFC,
                      pCutoff,
                      pAdjustment,
                      DESeq.args) {
  
  ## create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = meta,
                                        design = design)
  
  ## perform differential expression
  dds <- do.call(DESeq2::DESeq, c(list(dds), DESeq.args))
  
  ## reconstruct the user-defined contrast
  contrast <- strsplit(contrast, "-")[[1]]
  
  ## extract differential expression results for the appropriate contrast
  deRes <- DESeq2::results(dds,
                           contrast = c(group, contrast[1], contrast[2]),
                           pAdjustMethod = pAdjustment)
  deRes <- as.data.frame(deRes)
  
  ## add 'ID' column to differentially expressed results
  deRes$ID <- rownames(deRes)
  
  ## omit features with NA values
  deRes <- na.omit(deRes)
  
  ## extract normalized expression values
  normExpr <- DESeq2::counts(dds, normalized = TRUE)
  
  ## only retain valid features
  normExpr <- normExpr[rownames(normExpr) %in% deRes$ID, ]
  
  ## format differential expression table as required by MIRit
  deRes <- identifyColNames(deRes)
  
  ## select significant features
  sig <- rownames(deRes[abs(deRes$logFC) > logFC &
                          deRes$adj.P.Val < pCutoff, ])
  
  ## define the parameters used
  met <- paste("DESeq2 (p < ", pCutoff, ", correction: ",
               pAdjustment, ")", sep = "")
  
  ## create a list with DE results
  deList <- list(data = deRes,
                 significant = sig,
                 method = met,
                 deObject = dds,
                 normExpr = normExpr)
  
  ## return differential expression results
  return(deList)
  
}





## perform differential expression with limma-voom
voom.DE <- function(counts,
                    group,
                    contrast,
                    meta,
                    design,
                    logFC,
                    pCutoff,
                    pAdjustment,
                    useVoomWithQualityWeights,
                    filterByExpr.args,
                    calcNormFactors.args,
                    voom.args,
                    lmFit.args,
                    eBayes.args) {
  
  ## identify numerator level and reference level
  contrast <- strsplit(contrast, "-")[[1]]
  
  ## convert group variable to factor with specified reference level
  meta[, group] <- stats::relevel(factor(meta[, group]), ref = contrast[2])
  
  ## create edgeR object from counts
  features <- edgeR::DGEList(counts = counts,
                             group = meta[, group],
                             samples = meta)
  
  ## filter features based on expression
  keep <- do.call(edgeR::filterByExpr, c(list(features),
                                         filterByExpr.args))
  features <- features[keep, , keep.lib.sizes = FALSE]
  
  ## normalize feature counts (default with TMM)
  features <- do.call(edgeR::calcNormFactors,
                      c(list(features), calcNormFactors.args))
  
  ## extract normalized expression matrix
  normExpr <- edgeR::cpm(features, normalized.lib.sizes = TRUE, log = TRUE)
  
  ## design the model
  des <- stats::model.matrix(design, data = meta)
  
  ## apply voom transformation (with or without quality weights)
  if (useVoomWithQualityWeights == TRUE) {
    v <- do.call(limma::voomWithQualityWeights,
                 c(list(features, design = des), voom.args))
  } else {
    v <- do.call(limma::voom, c(list(features, design = des), voom.args))
  }
  
  ## fit a linear model for each gene
  fit <- do.call(limma::lmFit, c(list(v, design = des), lmFit.args))
  
  ## determine if the supplied model has intercept
  intercept <- attributes(terms(design))["intercept"] == 1
  
  ## identify the comparison for DE analysis
  if (intercept == FALSE) {
    
    ## build the contrast of interest
    contrast <- paste(group, contrast, sep = "")
    contrast <- paste(contrast, collapse = "-")
    
    ## fit the contrast
    con <- limma::makeContrasts(contrasts = contrast,
                                levels = des)
    
    ## set contrast parameter for DE
    comparison <- list(contrast = con)
    
  } else {
    
    ## set the coefficient name for the appropriate comparison
    cf <- contrast[1]
    cf <- paste(group, cf, sep = "")
    
    ## set contrast parameter for DE
    comparison <- list(coef = cf)
    
  }
  
  ## fit the contrast matrix
  fit <- do.call(limma::contrasts.fit, c(list(fit), comparison))
  
  ## perform empirical Bayes smoothing
  fit <- do.call(limma::eBayes, c(list(fit), eBayes.args))
  
  ## retrieve differentially expressed features
  compName <-ifelse(intercept == TRUE, comparison[[1]], contrast)
  deRes <- limma::topTable(fit,
                           coef = compName,
                           number = Inf,
                           adjust.method = pAdjustment)
  
  ## add 'ID' column to differentially expressed results
  deRes$ID <- rownames(deRes)
  
  ## format differential expression table as required by MIRit
  deRes <- identifyColNames(deRes)
  
  ## select significant features
  sig <- rownames(deRes[abs(deRes$logFC) > logFC &
                          deRes$adj.P.Val < pCutoff, ])
  
  ## define the parameters used
  met <- paste("limma-voom (p < ", pCutoff, ", correction: ",
               pAdjustment, ")", sep = "")
  
  ## create a list with DE results
  deList <- list(data = deRes,
                 significant = sig,
                 method = met,
                 deObject = features,
                 normExpr = normExpr)
  
  ## return differential expression results
  return(deList)
  
}





## perform differential expression with limma
limma.DE <- function(expr,
                     group,
                     contrast,
                     meta,
                     design,
                     logFC,
                     pCutoff,
                     pAdjustment,
                     useArrayWeights,
                     useWsva,
                     wsva.args,
                     arrayWeights.args,
                     lmFit.args,
                     eBayes.args) {
  
  ## identify numerator level and reference level
  contrast <- strsplit(contrast, "-")[[1]]
  
  ## convert group variable to factor with specified reference level
  meta[, group] <- stats::relevel(factor(meta[, group]), ref = contrast[2])
  
  ## design the linear model
  des <- stats::model.matrix(design, data = meta)
  
  ## correct for unknown sources of variability with WSVA
  if (useWsva == TRUE) {
    wcov <- do.call(limma::wsva, c(list(expr, design = des), wsva.args))
    des <- cbind(des, wcov)
  }
  
  ## fit a linear model for each gene with or without array quality weights
  if (useArrayWeights == TRUE) {
    arrayw <- do.call(limma::arrayWeights, c(list(expr), arrayWeights.args))
    fit <- do.call(limma::lmFit, c(list(expr, design = des, weights = arrayw),
                                   lmFit.args))
  } else {
    fit <- do.call(limma::lmFit, c(list(expr, design = des), lmFit.args))
  }
  
  ## determine if the supplied model has intercept
  intercept <- attributes(terms(design))["intercept"] == 1
  
  ## identify the comparison for DE analysis
  if (intercept == FALSE) {
    
    ## build the contrast of interest
    contrast <- paste(group, contrast, sep = "")
    contrast <- paste(contrast, collapse = "-")
    
    ## fit the contrast
    con <- limma::makeContrasts(contrasts = contrast,
                                levels = des)
    
    ## set contrast parameter for DE
    comparison <- list(contrast = con)
    
  } else {
    
    ## set the coefficient name for the appropriate comparison
    cf <- contrast[1]
    cf <- paste(group, cf, sep = "")
    
    ## set contrast parameter for DE
    comparison <- list(coef = cf)
    
  }
  
  ## fit the contrast matrix
  fit <- do.call(limma::contrasts.fit, c(list(fit), comparison))
  
  ## perform empirical Bayes smoothing
  fit <- do.call(limma::eBayes, c(list(fit), eBayes.args))
  
  ## retrieve differentially expressed features
  compName <-ifelse(intercept == TRUE, comparison[[1]], contrast)
  deRes <- limma::topTable(fit,
                           coef = compName,
                           number = Inf,
                           adjust.method = pAdjustment)
  
  ## add 'ID' column to differentially expressed results
  deRes$ID <- rownames(deRes)
  
  ## format differential expression table as required by MIRit
  deRes <- identifyColNames(deRes)
  
  ## select significant features
  sig <- rownames(deRes[abs(deRes$logFC) > logFC &
                          deRes$adj.P.Val < pCutoff, ])
  
  ## define the parameters used
  met <- paste("limma (p < ", pCutoff, ", correction: ",
               pAdjustment, ")", sep = "")
  
  ## create a list with DE results
  deList <- list(data = deRes,
                 significant = sig,
                 method = met,
                 deObject = fit,
                 normExpr = normExpr)
  
  ## return differential expression results
  return(deList)
  
}





#' Manually add differential expression results to a MirnaExperiment object
#' 
#' This function allows to add miRNA and gene differential expression results
#' to a [`MirnaExperiment`][MirnaExperiment-class] object. Instead of running
#' [performMirnaDE()] and [performGeneDE()] functions, this one allows to use
#' differential expression analyses carried out in other ways. This is
#' particularly useful in order to use the pipeline implemented in MIRit for
#' proteomic data and for expression data deriving from different technologies.
#' 
#' @details
#' The following paragraphs briefly explain the formats needed for mirnaDE,
#' geneDE, significantMirnas and significantGenes.
#' 
#' ## mirnaDE and geneDE
#'
#' `mirnaDE` and `geneDE` are two objects of class `data.frame` containing
#' the results of miRNA and gene differential expression analysis respectively.
#' These tables should contain the differential expression results for all
#' miRNAs/genes analyzed, not just for statistically significant species.
#' All `data.frame` objects can be used, as long as they have:
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
#' symbols in `significantGenes` must be present in the `ID` column of `geneDE`.
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
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
#' 
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing differential
#' expression results. To access these results, the user may run the
#' [mirnaDE()] and [geneDE()] functions for miRNAs and genes, respectively.
#' 
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # create samples metadata
#' meta <- data.frame("primary" = colnames(geneCounts),
#' "mirnaCol" = colnames(mirnaCounts), "geneCol" = colnames(geneCounts),
#' "disease" = c(rep("PTC", 8), rep("NTH", 8)), 
#' "patient" = c(rep(paste("Sample_", seq(8), sep = ""), 2)))
#' 
#' # perform miRNA DE with DESeq2 separately
#' dds_m <- DESeq2::DESeqDataSetFromMatrix(countData = obj[["microRNA"]],
#' colData = meta, design = ~ 0 + disease + patient)
#' dds_m <- DESeq2::DESeq(dds_m)
#' de_m <- as.data.frame(DESeq2::results(dds_m,
#' contrast = c("disease", "PTC", "NTH"), pAdjustMethod = "fdr"))
#' 
#' # perform gene DE with DESeq2 separately
#' dds_g <- DESeq2::DESeqDataSetFromMatrix(countData = obj[["genes"]],
#' colData = meta, design = ~ 0 + disease + patient)
#' dds_g <- DESeq2::DESeq(dds_g)
#' de_g <- as.data.frame(DESeq2::results(dds_g,
#' contrast = c("disease", "PTC", "NTH"), pAdjustMethod = "fdr"))
#' 
#' # prepare DE tables
#' de_m$ID <- rownames(de_m)
#' de_m <- na.omit(de_m)
#' de_g$ID <- rownames(de_g)
#' de_g <- na.omit(de_g)
#' 
#' # define significant features
#' sig_m <- de_m$ID[de_m$padj < 0.05, ]
#' sig_g <- de_g$ID[de_g$padj < 0.05, ]
#' 
#' # add DE results to MirnaExperiment object
#' obj <- addDifferentialExpression(obj, de_m, de_g, sig_m, sig_g)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#' 
#' @export
addDifferentialExpression <- function(mirnaObj,
                                      mirnaDE,
                                      geneDE,
                                      significantMirnas,
                                      significantGenes) {
  
  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
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
  
  ## check and set accepted column names in DE data.frames
  mirnaDE <- identifyColNames(mirnaDE, "miRNA")
  geneDE <- identifyColNames(geneDE, "gene")
  
  ## check that miRNA and gene names are the same across data
  if (!identical(sort(rownames(mirnaObj[["microRNA"]])), sort(mirnaDE$ID))) {
    stop(paste("Row names of 'mirnaObj[['microRNA']]' must match the miRNA",
               "identifiers present in 'mirnaDE'"), call. = FALSE)
  }
  if (!identical(sort(rownames(mirnaObj[["microRNA"]])), sort(geneDE$ID))) {
    stop(paste("Row names of 'mirnaObj[['genes']]' must match the gene",
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
  
  ## add miRNA and gene differential expression results to mirnaObj
  mirnaDE(mirnaObj) <- list(data = mirnaDE,
                            significant = significantMirnas,
                            method = "manually added",
                            deObject = NULL)
  geneDE(mirnaObj) <- list(data = geneDE,
                           significant = significantGenes,
                           method = "manually added",
                           deObject = NULL)
  
  ## return mirnaObj
  return(mirnaObj)
  
}

