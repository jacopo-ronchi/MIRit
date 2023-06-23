#' Integrate microRNA and gene expression
#'
#' This function allows to identify microRNAs that are significantly
#' associated/correlated with their targets. The principle is that, since the
#' biological role of miRNAs is mainly to negatively regulate gene expression
#' post-transcriptionally, the expression of a microRNA should be negatively
#' correlated with the expression of its targets. To test this assumption for
#' matched-sample data, this function performs a correlation analysis. On the
#' other hand, for unpaired data, it offers different one-sided association
#' tests to estimate if targets of down-regulated miRNAs are enriched in
#' up-regulated genes and vice versa. Additionally, for unpaired data, miRNA
#' effects on target gene expression can also be quantified through a fast
#' approximation to rotation gene-set testing ('fry' method). For correlation
#' analyses, the default behavior is to use Spearman's correlation analysis,
#' whereas for association tests the default option makes use of a one-sided
#' Fisher's exact test with Lancaster's mid-p correction. See the *details*
#' section for further information.
#'
#' @details
#' As already pointed out, if miRNA and gene expression data derive from the
#' same samples, a correlation analysis is used. For evaluating these
#' relationships, the default method used is Spearman's correlation
#' coefficient, as:
#' * it does not need normally distributed data;
#' * it does not assume linearity;
#' * it is much more resistant to outliers.
#'
#' However, the user can also decide to use other correlation methods,
#' such as Pearson's and Kendall's correlation. Nevertheless, for NGS data
#' it may happen that a certain number of ties is present in the expression
#' values. This can be handled by `spearman` method as it computes a
#' tie-corrected version of Spearman's coefficients. However, another
#' correlation method that is suitable to perform rank correlation on tied data
#' is the Kendall's tau-b method, usable with `kendall`.
#' 
#' Regarding correlation direction, since miRNAs mainly act as negative
#' regulators, only negatively correlated miRNA-target pairs are evaluated, and
#' statistical significance is calculated through a one-tailed t-test.
#'
#' Moreover, if gene expression data and miRNA expression data derive from
#' different samples (unpaired data), a correlation analysis can't be
#' performed. However, one-sided association tests can be applied in these
#' cases to evaluate if targets of down-regulated miRNAs are statistically
#' enriched in up-regulated genes, and, conversely, if targets of up-regulated
#' miRNAs are statistically enriched in down-regulated genes. In this case,
#' Fisher's exact test is the default option to assess the statistical
#' significance of this inverse association. However, if cases and controls
#' derive from repeated measurements taken in different times or from different
#' tissues, a one-sided McNemar's exact test is used instead. For both Fisher's
#' and McNemar's tests, the default behavior is to use Lancaster's mid-p
#' adjustment since it has been shown that it increases statistical power
#' while retaining Type I error rates.
#' 
#' Finally, for unpaired data, the effect of DE-miRNAs on the expression of
#' target genes can be estimated through rotation gene-set tests. In particular,
#' a fast approximation to rotation gene-set testing called `fry`, implemented
#' in the `limma` package, can be used to statistically quantify the influence
#' of miRNAs on the expression changes of theri target genes.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param test The statistical test to evaluate the association between miRNAs
#' and genes. It must be one of `auto` (default), to automatically determine
#' the appropriate statistical test; `correlation`, to perform a correlation
#' analysis; `association`, to perform a one-sided association test; `fry` to
#' perform the integrative analysis through rotation gene-set testing
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param corMethod The correlation method to be used for correlation analysis.
#' It must be one of: `spearman` (default), `pearson`, `kendall`. See the
#' *details* section for further information
#' @param corCutoff The minimum (negative) value of correlation coefficient to
#' consider meaningful a miRNA-target relationship. Default is `0.5`
#' @param midpAdjustment Logical, whether to use Lancaster's mid-p adjustment
#' method for association tests or not. Default is TRUE to compute mid-p-values
#' @param paired Logical, this parameter must be set to TRUE if samples from
#' the biological conditions considered for differential expression analysis
#' derive from repeated measurements. In this case, a one-sided McNemar's exact
#' test is used for estimating the association. Default is FALSE
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing integration
#' results. To access these results, the user can make use of the
#' [integration()] function. For additional details on how to
#' interpret the results of miRNA-gene integrative analysis, please see
#' [`MirnaExperiment`][MirnaExperiment-class].
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform integration analysis with default settings
#' obj <- mirnaIntegration(obj)
#'
#' # use Fisher's exact test  with FDR < 0.05 as significance threshold
#' obj <- mirnaIntegration(obj, test = "association",
#' pAdjustment = "fdr")
#'
#' # perform Kendall's correlation analysis with tau > 0.8 and p < 0.05
#' obj <- mirnaIntegration(obj, test = "correlation",
#' corMethod = "kendall", corCutoff = 0.8)
#' 
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies.” Nucleic Acids Research, 43(7), e47. \url{doi:10.1093/nar/gkv007}.
#' 
#' Di Wu and others, ROAST: rotation gene set tests for complex microarray
#' experiments, Bioinformatics, Volume 26, Issue 17, September 2010,
#' Pages 2176–2182, \url{https://doi.org/10.1093/bioinformatics/btq401}.
#' 
#' Routledge, R. D. (1994). Practicing Safe Statistics with the Mid-p. The
#' Canadian Journal of Statistics / La Revue Canadienne de Statistique, 22(1),
#' 103–110, \url{https://doi.org/10.2307/3315826}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
mirnaIntegration <- function(mirnaObj,
                             test = "auto",
                             pCutoff = 0.05,
                             pAdjustment = "fdr",
                             corMethod = "spearman",
                             corCutoff = 0.5,
                             midpAdjustment = TRUE,
                             paired = FALSE) {
  
  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("MiRNA differential expression results are not present in",
               "'mirnaObj'. Please, use 'performMirnaDE()' before using",
               "this function. See ?performMirnaDE"), call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (max(dim(mirnaTargets(mirnaObj))) == 0) {
    stop(paste("No targets are present within 'mirnaObj'!",
               "Before performing the integration analysis miRNA target",
               "genes must be retrieved with the 'getTargets()' function.",
               "See '?getTargets' for the details."), call. = FALSE)
  }
  if (!is.character(test) |
      length(test) != 1 |
      !test %in% c("auto", "correlation", "association", "fry")) {
    stop(paste("'test' must be one of:\n",
               "\t- 'auto', (default) to automatically",
               "choose the appropriate test;\n",
               "\t- 'correlation', to perform a correlation analysis;\n",
               "\t- 'association', to apply a one-sided association test;\n",
               "\t- 'fry', to apply rotation gene-set testing."),
         call. = FALSE)
  }
  if (test == "correlation" & pairedSamples(mirnaObj) == FALSE) {
    warning(paste("You can't perform a 'correlation' analysis with",
                  "unpaired samples! Setting 'test' to 'association' for",
                  "performing a one-sided association test instead..."),
            call. = FALSE)
    test <- "association"
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
  if (!is.character(corMethod) |
      length(corMethod) != 1 |
      !corMethod %in% c("pearson", "kendall", "spearman")) {
    stop(paste("'corMethod' must be one of: 'pearson', 'kendall' and",
               "'spearman' (default)."),
         call. = FALSE)
  }
  if (!is.numeric(corCutoff) |
      length(corCutoff) != 1 |
      corCutoff > 1 |
      corCutoff < 0) {
    stop("'corCutoff' must be a number between 0 and 1! (default is 0.5)",
         call. = FALSE)
  }
  
  ## use the appropriate test
  if (test == "auto" & pairedSamples(mirnaObj) == TRUE) {
    intMethod <- "correlation"
    message(paste("Since data derive from paired samples, a correlation test",
                  "will be used."))
  } else if (test == "auto" & pairedSamples(mirnaObj) == FALSE) {
    intMethod <- "association"
    message(paste("Since data derive from different samples, a",
                  "one-sided association test will be used."))
  } else {
    intMethod <- test
    message(paste("As specified by the user, a",
                  ifelse(test == "correlation", "correlation",
                         "one-sided association test"),
                  "will be used."))
  }
  
  ## call the right function
  if (intMethod == "correlation") {
    mirnaObj <- correlateMirnaTargets(mirnaObj,
                                      corMethod,
                                      corCutoff,
                                      pCutoff,
                                      pAdjustment)
  } else if (intMethod == "association") {
    mirnaObj <- associateMirnaTargets(mirnaObj,
                                      pCutoff,
                                      pAdjustment,
                                      midpAdjustment,
                                      paired)
  } else if (intMethod == "fry") {
    mirnaObj <- fryMirnaTargets(mirnaObj,
                                pCutoff,
                                pAdjustment)
  }
  
  ## return the object with integration slot
  return(mirnaObj)
  
}





## correlation analysis
correlateMirnaTargets <- function(mirnaObj,
                                  corMethod,
                                  corCutoff,
                                  pCutoff,
                                  pAdjustment) {
  
  ## extract miRNA and gene expression values
  mirnaExpr <- mirnaObj[["microRNA"]]
  geneExpr <- mirnaObj[["genes"]]
  
  ## check if samples are paired, otherwise exclude unpaired samples
  sMap <- MultiAssayExperiment::sampleMap(mirnaObj)
  mirnaSamples <- sMap$primary[sMap$assay == "microRNA"]
  geneSamples <- sMap$primary[sMap$assay == "genes"]
  
  if (!identical(mirnaSamples, geneSamples)) {
    
    ## determine common and uncommon samples
    common <- intersect(mirnaSamples, geneSamples)
    unpaired  <- setdiff(mirnaSamples, geneSamples)
    
    ## stop if common samples are less than 3
    if (length(common) < 3) {
      stop("To perform the correlation at least 3 paired samples are needed.",
           call. = FALSE)
    }
    
    ## remove samples without measurments of both miRNAs and genes
    if (length(unpaired) > 0) {
      
      mirnaExpr <- mirnaExpr[, sMap$colname[sMap$assay == "microRNA" &
                                              sMap$primary %in% common]]
      geneExpr <- geneExpr[, sMap$colname[sMap$assay == "genes" &
                                            sMap$primary %in% common]]
      
      warning(paste("Some samples don't have expression values for both",
                    "miRNAs and genes, thus they have been excluded from the",
                    "correlation analysis. Removed samples are:",
                    paste(unpaired, collapse = ", ")))
    }
    
    ## order the columns of expression matrices in the same way
    mirnaMap = sMap[sMap$assay == "microRNA", ]
    geneMap = sMap[sMap$assay == "genes", ]
    mirnaOrder <- mirnaMap$primary[order(match(mirnaMap$colname,
                                               colnames(mirnaExpr)))]
    geneExpr <- geneExpr[, geneMap$colname[order(match(geneMap$primary,
                                                       mirnaOrder))]]
    
  }
  
  ## retrieve differentially expressed genes and miRNAs
  dem <- mirnaDE(mirnaObj)
  deg <- geneDE(mirnaObj)
  
  ## retrieve targets of DE-miRNAs from the object
  targetsTable <- mirnaTargets(mirnaObj)
  
  ## select differentially expressed miRNA targets
  targetsTable <- targetsTable[targetsTable$Gene.Symbol %in% deg$ID, ]
  
  ## restrict to target genes present in the assay
  targetsTable <- targetsTable[targetsTable$Gene.Symbol
                               %in% rownames(geneExpr), ]
  
  ## extract the expression values of miRNA targets
  targetExpr <- geneExpr[rownames(geneExpr) %in% targetsTable$Gene.Symbol, ]
  
  ## compute the correlation between each pair of DE-miRNA - target
  usedCoef <- gsub("(^)([[:alpha:]])", "\\1\\U\\2", corMethod, perl = TRUE)
  message(paste("Performing ", usedCoef, "'s correlation analysis...",
                sep = ""))
  correlation <- mapply(function(mirna, gene) {
    
    ## extract the expression values of the microRNA and the target
    mirnaInt <- as.numeric(mirnaExpr[mirna, ])
    geneInt <- as.numeric(targetExpr[gene, ])
    
    ## perform the correlation analysis
    corPair <- stats::cor.test(mirnaInt,
                               geneInt,
                               method = corMethod,
                               alternative = "less",
                               exact = FALSE)
    
    ## report the results of the correlation analysis
    fold <- ifelse(dem$logFC[dem$ID == mirna] > 0,
                   "upregulated",
                   "downregulated")
    
    pair <- c(mirna,
              gene,
              fold,
              corPair$estimate,
              corPair$p.value)
    pair
    
  }, targetsTable$MicroRNA, targetsTable$Gene.Symbol)
  
  ## convert correlation output to a data.frame object
  corRes <- as.data.frame(t(correlation))
  colnames(corRes) <-c("microRNA", "Target", "microRNA.Direction",
                       "Corr.Coefficient", "Corr.P.Value")
  corRes$Corr.Coefficient <- as.numeric(corRes$Corr.Coefficient)
  corRes$Corr.P.Value <- as.numeric(corRes$Corr.P.Value)
  
  ## correct correlation p values for multiple testing
  pAdj <- stats::p.adjust(corRes$Corr.P.Value, method = pAdjustment)
  corRes$Corr.Adjusted.P.Val <- pAdj
  
  ## select statistically significant associations
  corRes <- corRes[corRes$Corr.Adjusted.P.Val <= pCutoff &
                     abs(corRes$Corr.Coefficient) >= corCutoff, ]
  
  ## report the results of the correlation analysis
  if (nrow(corRes) >= 1) {
    message(paste("A statistically significant correlation between",
                  nrow(corRes), "miRNA-target pairs was found!"))
  } else {
    message(paste("No statistically significant correlation between",
                  "DE-miRNAs and targets was found."))
  }
  
  ## return the object with the results of the correlation analysis
  resList <- list(data = corRes,
                  method = paste(usedCoef, "'s correlation analysis", sep = ""),
                  pCutoff = pCutoff,
                  pAdjustment = pAdjustment)
  integration(mirnaObj) <- resList
  return(mirnaObj)
  
}





## one-sided association tests for integrative analysis
associateMirnaTargets <- function(mirnaObj,
                                  pCutoff,
                                  pAdjustment,
                                  midpAdjustment,
                                  paired) {
  
  ## obtain differentially expressed miRNAs and genes
  dem <- mirnaDE(mirnaObj)
  deg <- geneDE(mirnaObj)
  
  ## divide DEGs in up- and down-regulated
  upDeg <- deg$ID[deg$logFC > 0]
  downDeg <- deg$ID[deg$logFC < 0]
  
  ## extract target genes and total genes
  targetsTable <- mirnaTargets(mirnaObj)
  targets <- targetsTable$Gene.Symbol
  targets <- unique(targets)
  totGenes <- rownames(mirnaObj[["genes"]])
  totGenes <- unique(totGenes)
  
  ## define the association test used
  if (paired == FALSE & midpAdjustment == FALSE) {
    meth <- "One-sided Fisher's exact test"
  } else if (paired == FALSE & midpAdjustment == TRUE) {
    meth <- paste("One-sided Fisher's exact test",
                  "with Lancaster's mid-p correction")
  } else if (paired == TRUE & midpAdjustment == FALSE) {
    meth <- "One-sided McNemar's exact test"
  } else if (paired == TRUE & midpAdjustment == TRUE) {
    meth <- paste("One-sided McNemar's exact test",
                  "with Lancaster's mid-p correction")
  }
  
  ## compute the inverse association between each miRNA and its targets
  message(paste("Performing ", meth, "...", sep = ""))
  association <- lapply(dem$ID, function(x) {
    
    ## differentiate between up-miRNA - down-genes and down-miRNA - up-genes
    fold <- ifelse(dem$logFC[dem$ID == x] > 0, "Up", "Down")
    
    if (fold == "Up") {
      degFold <- downDeg
      geneFold <- "Down"
    } else if (fold == "Down") {
      degFold <- upDeg
      geneFold <- "Up"
    }
    
    ## build the contingency table
    degs <- degFold
    nonDegs <- setdiff(totGenes, degs)
    nonTarg <- setdiff(
      totGenes,
      targetsTable$Gene.Symbol[targetsTable$MicroRNA == x]
    )
    targ <- intersect(
      totGenes,
      targetsTable$Gene.Symbol[targetsTable$MicroRNA == x]
    )
    
    nonDegsNonTarg <- length(intersect(nonDegs, nonTarg))
    degsNonTarg <- length(intersect(degs, nonTarg))
    nonDegsTarg <- length(intersect(nonDegs, targ))
    degsTarg <- length(intersect(degs, targ))
    
    contingencyTable <- data.frame(
      "NonTargetGenes" = c(nonDegsNonTarg, degsNonTarg),
      "TargetGenes" = c(nonDegsTarg, degsTarg),
      row.names = c("nonDEGs", "DEGs"),
      stringsAsFactors = FALSE
    )
    
    ## perform the appropriate statistical test
    pval <- association.helper(contingencyTable,
                               midpAdjustment,
                               paired)
    
    ## report results to a list object
    intRes <- c(x,
                fold,
                geneFold,
                degsTarg,
                degsTarg + nonDegsTarg,
                pval,
                paste(intersect(degs, targ), collapse = "/"))
    intRes
    
  })
  
  ## convert list object into a data.frame
  association <- as.data.frame(do.call(rbind, association))
  colnames(association) <- c("microRNA",
                             "mirna.direction",
                             "gene.direction",
                             "DE",
                             "targets",
                             "P.Val",
                             "DE.targets")
  
  ## correct p-values
  pAdj <- stats::p.adjust(association$P.Val, method = pAdjustment)
  association$adj.P.Val <- pAdj
  association <- association[, c(1, 2, 3, 4, 5, 6, 8, 7)]
  association <- association[order(association$adj.P.Val), ]
  
  ## print association results
  association <- association[association$adj.P.Val <= pCutoff, ]
  if (nrow(association) >= 1) {
    message(paste("A statistically significant association between",
                  nrow(association), "DE-miRNAs and",
                  length(unique(unlist(strsplit(association$DE.targets, "/")))),
                  "genes was found!"))
  } else {
    message(paste("No statistically significant associations between",
                  "DE-miRNAs and DEGs were found."))
  }
  
  ## return the object containing association results
  resList <- list(data = association,
                  method = meth,
                  pCutoff = pCutoff,
                  pAdjustment = pAdjustment)
  integration(mirnaObj) <- resList
  return(mirnaObj)
  
}





## helper function for choosing the association test
association.helper <- function(contingencyTable, midpAdjustment, paired) {
  
  ## perform the appropriate statistical test
  if (paired == FALSE) {
    ## perform one-sided Fisher's test with/without Lancaster's correction
    pval <- fisher.midp(contingencyTable, midpAdjustment)
  } else {
    ## perform one-sided McNemar's test with/without Lancaster's correction
    pval <- mcnemar.midp(contingencyTable, midpAdjustment)
  }
  
  ## return the computed p-value
  return(pval)
  
}





## compute p-value for one-sided Fisher's exact test with or without
## Lancaster's mid-p correction (alternative is greater)
fisher.midp <- function(mat, midpAdjustment) {
  m <- sum(mat[, 1L])
  n <- sum(mat[, 2L])
  k <- sum(mat[1L, ])
  mat <- mat[1L, 1L]
  if (midpAdjustment == TRUE) {
    midp <- stats::phyper(mat - 1, m, n, k,
                          lower.tail = FALSE) - 0.5 * stats::dhyper(mat, m,
                                                                    n, k)
  } else {
    midp <- stats::phyper(mat - 1, m, n, k, lower.tail = FALSE)
  }
  return(midp)
}





## compute p-value for one-sided McNemar's exact test with or without
## Lancaster's mid-p correction (alternative is greater)
mcnemar.midp <- function(mat, midpAdjustment) {
  x <- mat[1, 2]
  n <- mat[1, 2] + mat[2, 1]
  if (n != 0){
    if (midpAdjustment == TRUE) {
      midp <- stats::pbinom(x-1, n, 0.5, lower.tail = FALSE) - 
        0.5 * stats::dbinom(x, n, 0.5)
    } else {
      midp <- stats::pbinom(x-1, n, 0.5, lower.tail = FALSE)
    }
  } else {
    midp <- 1
  }
  return(midp)
}





## rotation gene-set test for integrative analysis
fryMirnaTargets <- function(mirnaObj,
                            pCutoff,
                            pAdjustment) {
  
  ## extract gene differential expression results
  de <- geneDE(mirnaObj, param = TRUE)
  
  ## access sample metadata
  meta <- MultiAssayExperiment::colData(mirnaObj)
  meta <- meta[!is.na(meta$geneCol), ]
  
  ## build miRNA-target sets
  tg <- mirnaTargets(mirnaObj)
  tgList <- split(tg$Gene.Symbol, tg$MicroRNA)
  
  ## determine the appropriate expression matrix and the experimental design
  if (de$method == "limma") {
    expr <- mirnaObj[["genes"]]
    des <- stats::model.matrix(de$design, data = meta)
  } else if (de$method == "edgeR" |
             de$method == "limma-voom") {
    expr <- geneDE(mirnaObj, returnObject = TRUE)
    des <- expr$design
  } else if (de$method == "DESeq2") {
    message("Applying 'limma-voom' pipeline before using 'fry' method...")
    counts <- MultiAssayExperiment::metadata(mirnaObj)[["oldCounts"]][["genes"]]
    features <- edgeR::DGEList(counts = counts,
                               group = meta[, de$group],
                               samples = meta)
    keep <- edgeR::filterByExpr(features)
    features <- features[keep, , keep.lib.sizes = FALSE]
    features <- edgeR::calcNormFactors(features)
    des <- stats::model.matrix(de$design, data = meta)
    expr <- limma::voom(features, design = des)
  }
  
  ## set up the contrast
  contrast <- strsplit(de$contrast, "-")[[1]]
  
  ## determine if the model has intercept
  intercept <- attributes(terms(de$design))["intercept"] == 1
  
  ## identify the comparison for DE analysis
  if (intercept == FALSE) {
    
    ## build the contrast of interest
    contrast <- paste(de$group, contrast, sep = "")
    contrast <- paste(contrast, collapse = "-")
    
    ## create contrast matrix
    con <- limma::makeContrasts(contrasts = contrast,
                                levels = des)
    
  } else {
    
    ## set the coefficient name for the appropriate comparison
    con <- contrast[1]
    con <- paste(de$group, con, sep = "")
    
  }
  
  ## perform integration through rotation gene set enrichment analysis
  message("Performing miRNA-gene integrative analysis using 'fry' method...")
  if (de$method == "edgeR") {
    
    ## perform miRNA-target integration through 'fry'
    rs <- edgeR::fry.DGEList(y = expr,
                             index = tgList,
                             design = des,
                             contrast = con,
                             adjust.method = pAdjustment)
    
  } else {
    
    ## perform miRNA-target integration through 'fry'
    rs <- limma::fry(y = expr,
                     index = tgList,
                     design = des,
                     contrast = con,
                     adjust.method = pAdjustment)
    
  }
  
  ## retain effects in the right direction
  dem <- mirnaDE(mirnaObj)
  upDem <- dem$ID[dem$logFC > 0]
  downDem <- dem$ID[dem$logFC < 0]
  rs <- rs[(rownames(rs) %in% upDem & rs$Direction == "Down") |
             (rownames(rs) %in% downDem & rs$Direction == "Up"), ]
  
  ## maintain interactions under the specified cutoff
  res <- rs[rs$FDR <= pCutoff, ]
  
  ## create resulting data.frame
  res$microRNA <- rownames(res)
  res$mirna.direction <- "Down"
  res$mirna.direction[res$microRNA %in% upDem] <- "Up"
  deg <- geneDE(mirnaObj)
  deTarg <- mapply(function(mir, fold) {
    if (fold == "Up") {
      degFold <- deg$ID[deg$logFC < 0]
    } else if (fold == "Down") {
      degFold <- deg$ID[deg$logFC > 0]
    }
    mirDe <- intersect(tgList[[mir]], degFold)
    c(paste(mirDe, collapse = "/"), length(mirDe))
  }, res$microRNA, res$mirna.direction)
  res$DE_targets <- deTarg[1, ]
  res$DE <- deTarg[2, ]
  res <- res[, c(7, 8, 2, 10, 1, 3, 4, 9)]
  colnames(res) <- c("microRNA", "mirna.direction", "gene.direction",
                     "DE", "targets", "P.Val", "adj.P.Val", "DE.targets")
  
  ## print integration results
  if (nrow(res) >= 1) {
    message(paste("A statistically significant association between",
                  nrow(res), "DE-miRNAs and",
                  length(unique(unlist(strsplit(res$DE.targets, "/")))),
                  "genes was found!"))
  } else {
    message(paste("No statistically significant associations between",
                  "DE-miRNAs and genes were found."))
  }
  
  ## return the object containing association results
  resList <- list(data = res,
                  method = "Rotation Gene-Set Test (FRY)",
                  pCutoff = pCutoff,
                  pAdjustment = pAdjustment)
  integration(mirnaObj) <- resList
  return(mirnaObj)
  
}

