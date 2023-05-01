#' Integrate microRNA and gene expression
#'
#' This function allows to test which microRNAs are significantly
#' associated/correlated with their targets. The principle is that, since the
#' biological role of miRNAs is mainly to negatively regulate gene expression
#' post-transcriptionally, the expression of a microRNA should be negatively
#' correlated with the expression of its targets. To test this assumption, this
#' function implements a correlation analysis, if miRNA and gene expression
#' data derive from the same samples (paired data), whereas, with unpaired data,
#' it uses a one-sided Fisher's exact to estimate if targets of down-regulated
#' miRNAs are enriched in up-regulated genes and vice versa. See the *details*
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
#' Finally, if gene expression data and miRNA expression data derive from
#' different samples (unpaired data), a correlation analysis can't be
#' performed. However, a one-sided Fisher's exact test can be applied in these
#' cases to evaluate if targets of down-regulated miRNAs are statistically
#' enriched in up-regulated genes, and, conversely, if targets of up-regulated
#' miRNAs are statistically enriched in down-regulated genes.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param test The statistical test to evaluate the association between miRNAs
#' and genes. It must be one of `auto` (default), to automatically determine
#' the appropriate statistical test; `correlation`, to perform a correlation
#' analysis; `fisher`, to perform a one-sided Fisher's exact test
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param onlySignificant Logical, whether to report only statistically
#' significant associations (default is `TRUE`). This parameter is considered
#' only if a one-sided Fisher's exact test is applied
#' @param corMethod The correlation method to be used for correlation analysis.
#' It must be one of: `spearman` (default), `pearson`, `kendall`. See the
#' 'details' section for further information
#' @param corCutoff The minimum (negative) value of correlation coefficient to
#' consider meaningful a miRNA-target relationship. Default is `0.5`
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing integration
#' results. To access these results, the user can make use of the
#' [mirnaTargetsIntegration()] function, which returns a `data.frame` object
#' with the details of miRNA-targets relationships. This `data.frame` differs
#' on the basis of the integration strategy used. For the one-sided Fisher's
#' exact test integration, it has seven columns:
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
#' * `DE_targets`: contains the list of differentially expressed targtes whose
#' expression is negatively associated with miRNA expression.
#' Instead, when a correlation analysis is performed, `mirnaTargetsIntegration`
#' has seven columns:
#' * `microRNA`: the miRNA ID;
#' * `target`: the correlated target gene;
#' * `microRNA.Direction`: the fold change direction of the DE-miRNA;
#' * `Corr.Coefficient`: the value of the correlation coefficient used;
#' * `Corr.P.Value`: the p-value resulting from the correlation analysis;
#' * `Corr.Adjusted.P.Val`: contains the correlation p-values corrected for
#' multiple testing.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform integration analysis with default settings
#' obj <- integrateMirnaTargets(obj)
#'
#' # use the Fisher's exact test  with FDR < 0.05 as significance threshold
#' obj <- integrateMirnaTargets(obj, test = "fisher", pAdjustment = "fdr")
#'
#' # perform Kendall's correlation analysis with tau > 0.8 and p < 0.05
#' obj <- integrateMirnaTargets(obj, test = "correlation",
#' corMethod = "kendall", corCutoff = 0.8)
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
integrateMirnaTargets <- function(mirnaObj,
                                  test = "auto",
                                  pCutoff = 0.05,
                                  pAdjustment = "fdr",
                                  onlySignificant = TRUE,
                                  corMethod = "spearman",
                                  corCutoff = 0.5) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (max(dim(mirnaTargets(mirnaObj))) == 0) {
    stop(paste("No targets are present within 'mirnaObj'!",
               "Before performing the integration analysis miRNA target",
               "genes must be retrieved with the 'getTargets()' function.",
               "See '?getTargets' for the details."), call. = FALSE)
  }
  if (!is.character(test) |
      length(test) != 1 |
      !test %in% c("auto", "correlation", "fisher")) {
    stop(paste("'test' must be one of:\n",
               "\t- 'auto', (default) to automatically",
               "choose the appropriate test;\n",
               "\t- 'correlation', to perform a correlation analysis;\n",
               "\t- 'fisher', to apply a one-sided Fisher's exact test."),
         call. = FALSE)
  }
  if (test == "correlation" & pairedSamples(mirnaObj) == FALSE) {
    warning(paste("You can't perform a 'correlation' analysis with",
                  "unpaired samples! Setting 'test' to 'fisher' for",
                  "performing a one-sided Fisher's exact test instead..."),
            call. = FALSE)
    test <- "fisher"
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
    stop(paste("'pAdjustment' must be  one of: 'none' (default), 'fdr',",
               "'BH' (same as 'fdr'), 'bonferroni', 'BY', 'hochberg',",
               "'holm', 'hommel'"),
         call. = FALSE)
  }
  if (!is.logical(onlySignificant) |
      length(onlySignificant) != 1) {
    stop("'onlySignificant' must be logical (TRUE/FALSE)!", call. = FALSE)
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
    intMethod <- "fisher"
    message(paste("Since data derive from different samples, a",
                  "one-sided Fisher's exact test will be used."))
  } else {
    intMethod <- test
    message(paste("As specified by the user, a",
                  ifelse(test == "correlation", "correlation",
                         "one-sided Fisher's exact test"),
                  "will be used."))
  }

  ## call the right function
  if (intMethod == "correlation") {
    mirnaObj <- correlateMirnaTargets(mirnaObj,
                                      corMethod,
                                      corCutoff,
                                      pCutoff,
                                      pAdjustment)
  } else if (intMethod == "fisher") {
    mirnaObj <- oneSidedFisher(mirnaObj, pCutoff, pAdjustment, onlySignificant)
  }

  ## return the object with integration slot
  return(mirnaObj)

}





## one-sided Fisher's exact test for association analysis
oneSidedFisher <- function(mirnaObj, pCutoff, pAdjustment, onlySignificant) {

  ## obtain differentially expressed miRNAs and genes
  dem <- mirnaDE(mirnaObj)
  deg <- geneDE(mirnaObj)

  ## divide DEGs in up- and down-regulated
  upDeg <- deg$ID[deg$logFC > 0]
  downDeg <- deg$ID[deg$logFC < 0]

  ## extract target genes and total genes
  targetsTable <- mirnaTargets(mirnaObj)
  targets <- targetsTable$target_symbol
  targets <- unique(targets)
  totGenes <- rownames(mirnaObj[["genes"]])
  totGenes <- unique(totGenes)

  ## compute the inverse association between each miRNA and its targets
  association <- lapply(dem$ID, function(x) {

    ## differentiate between up-miRNA - down-genes and down-miRNA - up-genes
    fold <- ifelse(dem$logFC[dem$ID == x] > 0, "upregulated", "downregulated")

    if (fold == "upregulated") {
      degFold <- downDeg
    } else if (fold == "downregulated") {
      degFold <- upDeg
    }

    ## build the contingency table
    degs <- degFold
    nonDegs <- setdiff(totGenes, degs)
    nonTarg <- setdiff(
      totGenes,
      targetsTable$target_symbol[targetsTable$mature_mirna_id == x]
      )
    targ <- intersect(
      totGenes,
      targetsTable$target_symbol[targetsTable$mature_mirna_id == x]
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

    ## perform one-sided Fisher's exact test
    fisherRes <- stats::fisher.test(contingencyTable, alternative = "greater")
    fisherPval <- fisherRes$p.value

    ## report results to a list object
    intRes <- c(x,
                fold,
                degsTarg,
                nonDegsTarg,
                fisherPval,
                paste(intersect(degs, targ), collapse = "/"))
    intRes

  })

  ## convert list object into a data.frame
  association <- as.data.frame(do.call(rbind, association))
  colnames(association) <- c("microRNA",
                             "direction",
                             "n_DE_targets",
                             "n_NON_DE_targets",
                             "Fisher.P.Val",
                             "DE_targets")

  ## correct p values
  adjFisher <- stats::p.adjust(association$Fisher.P.Val, method = pAdjustment)
  association$Fisher.Adjusted.P.Val <- adjFisher
  association <- association[, c(1, 2, 3, 4, 5, 7, 6)]
  association <- association[order(association$Fisher.Adjusted.P.Val), ]

  ## print association results
  if (onlySignificant == TRUE) {
    association <- association[association$Fisher.Adjusted.P.Val < pCutoff, ]
    if (nrow(association) >= 1) {
      message(paste("A statistically significant association between",
                    nrow(association), "DE-miRNAs and",
                    length(unique(unlist(stringr::str_split(
                      paste(association$DE_targets, collapse = "/"),
                      pattern = "/")))), "genes was found!"))
    } else {
      message(paste("No statistically significant associations between",
                    "DE-miRNAs and DEGs were found."))
    }
  } else {
    message("To access the results, use 'mirnaTargetsIntegration()'")
  }

  ## return the object containing association results
  mirnaTargetsIntegration(mirnaObj) <- association
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
  mirnaSamples <- colnames(mirnaExpr)
  geneSamples <- colnames(geneExpr)
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
      mirnaExpr <- mirnaExpr[, colnames(mirnaExpr) %in% common]
      geneExpr <- geneExpr[, colnames(geneExpr) %in% common]
      warning(paste("Some samples don't have expression values for both",
                    "miRNAs and genes, thus they have been excluded from the",
                    "correlation analysis. Removed samples are:",
                    paste(unpaired, collapse = ", ")))
    }

    ## order the columns of geneExpr according to those of mirnaExpr
    geneExpr <- geneExpr[, colnames(mirnaExpr)]

  }

  ## retrieve differentially expressed genes and miRNAs
  dem <- mirnaDE(mirnaObj)
  deg <- geneDE(mirnaObj)

  ## retrieve targets of DE-miRNAs from the object
  targetsTable <- mirnaTargets(mirnaObj)

  ## select differentially expressed miRNA targets
  targetsTable <- targetsTable[targetsTable$target_symbol %in% deg$ID, ]

  ## restrict to target genes present in the assay
  targetsTable <- targetsTable[targetsTable$target_symbol
                               %in% rownames(geneExpr), ]

  ## extract the expression values of miRNA targets
  targetExpr <- geneExpr[rownames(geneExpr) %in% targetsTable$target_symbol, ]

  ## compute the correlation between each pair of DE-miRNA - target
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

  }, targetsTable$mature_mirna_id, targetsTable$target_symbol)

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
  corRes <- corRes[corRes$Corr.Adjusted.P.Val < pCutoff &
                     abs(corRes$Corr.Coefficient) > corCutoff, ]

  ## report the results of the correlation analysis
  if (nrow(corRes) >= 1) {
    message(paste("A statistically significant correlation between",
                  nrow(corRes), "miRNA-target pairs was found!"))
  } else {
    message(paste("No statistically significant correlation between",
                  "DE-miRNAs and targets was found."))
  }

  ## return the object with the results of the correlation analysis
  mirnaTargetsIntegration(mirnaObj) <- corRes
  return(mirnaObj)

}

