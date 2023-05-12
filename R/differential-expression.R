performDE <- function(mirnaObj,
                      group,
                      contrast,
                      design,
                      assay = "miRNAs",
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
               "(e.g. 'PTC-NTH'). For details, see ?performDE"),
         call. = FALSE)
  }
  if (!rlang::is_formula(design) |
      !all(labels(stats::terms(design)) %in%
           colnames(MultiAssayExperiment::colData(mirnaObj)) &
           !labels(stats::terms(design)) %in%
           c("primary", "mirnaCol", "geneCol"))) {
    stop(paste("'design' must be an R formula that specifies the variables",
               "in colData that will be used to model expression. For",
               "details, see ?performDE"),
         call. = FALSE)
  }
  if (!is.character(assay) |
      length(assay) != 1 |
      !assay %in% c("miRNAs", "genes")) {
    stop(paste("'assay' must be either 'miRNAs' or 'genes'.",
               "For additional details, see ?performDE"),
         call. = FALSE)
  }
  if (!is.character(method) |
      length(method) != 1 |
      !method %in% c("limma", "edgeR", "DESeq2", "voom")) {
    stop(paste("'method' must be  one of: 'limma', 'edgeR' (default),",
               "'DESeq2', 'voom'. For additional details, see ?performDE"),
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
               "functions must be passed as lists! See ?performDE"),
         call. = FALSE)
  }
  
  ## define assay name and any previous DE object
  if (assay == "miRNAs") {
    assayName <- "microRNA"
    assayFunc <- "'mirnaDE()'"
    prevDe <- mirnaDE(mirnaObj, returnObject = TRUE)
  } else if (assay == "genes") {
    assayName <- "genes"
    assayFunc <- "'geneDE()'"
    prevDe <- geneDE(mirnaObj, returnObject = TRUE)
  }
  
  ## check if differential expression has already been carried out
  if (!is.null(prevDe)) {
    
    ## remove previous results and set expression back to counts
    if (assay == "miRNAs") {
      mirnaDE(mirnaObj) <- list(data = data.frame(),
                                significant = character(),
                                method = character(),
                                deObject = NULL)
      mirnaObj[[assayName]] <-
        MultiAssayExperiment::metadata(mirnaObj)[["oldCounts"]][[assayName]]
    } else if (assay == "genes") {
      geneDE(mirnaObj) <- list(data = data.frame(),
                               significant = character(),
                               method = character(),
                               deObject = NULL)
      mirnaObj[[assayName]] <-
        MultiAssayExperiment::metadata(mirnaObj)[["oldCounts"]][[assayName]]
    }
    
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
  if (assay == "miRNAs") {
    meta <- samplesMetadata[which(!is.na(samplesMetadata$mirnaCol)), ]
  } else if (assay == "genes") {
    meta <- samplesMetadata[which(!is.na(samplesMetadata$geneCol)), ]
  }
  
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
    method <- "limma"
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
  
  ## reconstruct the user-defined contrast
  contrast <- strsplit(contrast, "-")[[1]]
  contrast <- paste(group, contrast, sep = "")
  contrast <- paste(contrast, collapse = "-")
  
  ## set the contrast of interest
  con <- limma::makeContrasts(contrasts = contrast,
                              levels = des)
  
  ## estimate dispersion and fit negative binomial distribution
  features <- do.call(edgeR::estimateDisp,
                      c(list(features, des), estimateDisp.args))
  fit <- do.call(edgeR::glmQLFit, c(list(features, des), glmQLFit.args))
  
  ## perform differential expression
  de <- do.call(edgeR::glmQLFTest,
                c(list(fit, contrast = con), glmQLFTest.args))
  
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
  
  ## reconstruct the user-defined contrast
  contrast <- strsplit(contrast, "-")[[1]]
  contrast <- paste(group, contrast, sep = "")
  contrast <- paste(contrast, collapse = "-")
  
  ## set the contrast of interest
  con <- limma::makeContrasts(contrasts = contrast,
                              levels = des)
  
  ## fit the contrast matrix
  fit <- limma::contrasts.fit(fit, contrasts = con)
  
  ## perform empirical Bayes smoothing
  fit <- do.call(limma::eBayes, c(list(fit), eBayes.args))
  
  ## retrieve differentially expressed features
  deRes <- limma::topTable(fit,
                           coef = contrast,
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
  
  ## reconstruct the user-defined contrast
  contrast <- strsplit(contrast, "-")[[1]]
  contrast <- paste(group, contrast, sep = "")
  contrast <- paste(contrast, collapse = "-")
  
  ## set the contrast of interest
  con <- limma::makeContrasts(contrasts = contrast,
                              levels = des)
  
  ## fit the contrast matrix
  fit <- limma::contrasts.fit(fit, contrasts = con)
  
  ## perform empirical Bayes smoothing
  fit <- do.call(limma::eBayes, c(list(fit), eBayes.args))
  
  ## retrieve differentially expressed features
  deRes <- limma::topTable(fit,
                           coef = contrast,
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

