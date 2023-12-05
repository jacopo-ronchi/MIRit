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
#' Boschloo's exact test. See the *details* section for further information.
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
#' Please notice that if strong batch effects are noticed in expression data,
#' it is recommended to remove them through the [batchCorrection()] function
#' implemented in MIRit.
#'
#' Moreover, if gene expression data and miRNA expression data derive from
#' different samples (unpaired data), a correlation analysis can't be
#' performed. However, one-sided association tests can be applied in these
#' cases to evaluate if targets of down-regulated miRNAs are statistically
#' enriched in up-regulated genes, and, conversely, if targets of up-regulated
#' miRNAs are statistically enriched in down-regulated genes. In this case,
#' Fisher's exact test can be used to assess the statistical significance of
#' this inverse association. Moreover, Lancaster's mid-p adjustment can be
#' applied since it has been shown that it increases statistical power
#' while retaining Type I error rates. However, Fisher's exact test is a
#' conditional test that requires the sum of both rows and columns of a
#' contingency table to be fixed. Notably, this is not true for genomic data
#' because it is likely that different datasets may lead to a different number
#' of DEGs. Therefore, the default behavior in MIRit is to use a variant of
#' Barnard's exact test, named Boschloo's exact test, that is suitable when
#' group sizes of contingency tables are variable. Moreover, it is possible to
#' demonstrate that Boschloo's test is uniformly more powerful compared to
#' Fisher's exact test.
#'
#' Finally, for unpaired data, the effect of DE-miRNAs on the expression of
#' target genes can be estimated through rotation gene-set tests. In particular,
#' a fast approximation to rotation gene-set testing called `fry`, implemented
#' in the `limma` package, can be used to statistically quantify the influence
#' of miRNAs on the expression changes of their target genes.
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
#' @param associationMethod The statistical test used for evaluating the
#' association between miRNAs and their targets for unpaired data. It must be
#' one of `boschloo` (default), to perform a one-sided Boschloo's exact test;
#' `fisher-midp`, to compute a one-sided Fisher's exact test with Lancaster's
#' mid-p correction; `fisher`, to perform a one-sided Fisher's exact test
#' @param nuisanceParam The number of nuisance parameter values considered for
#' p-value calculation in `boschloo` method. The higher this value, the better
#' the p-value estimation accuracy. Default is 100
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
#' # perform Kendall's correlation analysis with tau > 0.8 and p < 0.05
#' obj <- mirnaIntegration(obj,
#'     test = "correlation",
#'     corMethod = "kendall", corCutoff = 0.8
#' )
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
#' Boschloo R.D. (1970). "Raised Conditional Level of Significance for the
#' 2x2-table when Testing the Equality of Two Probabilities".
#' Statistica Neerlandica. 24: 1–35.
#' \url{doi:10.1111/j.1467-9574.1970.tb00104.x}.
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
                             associationMethod = "boschloo",
                             nuisanceParam = 100) {
    ## check inputs
    if (!is(mirnaObj, "MirnaExperiment")) {
        stop("'mirnaObj' should be of class MirnaExperiment! ",
             "See ?MirnaExperiment",
             call. = FALSE
        )
    }
    if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0) {
        stop("MiRNA differential expression results are not present in ",
             "'mirnaObj'. Please, use 'performMirnaDE()' before using ",
             "this function. See ?performMirnaDE",
             call. = FALSE
        )
    }
    if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
        stop("Gene differential expression results are not present in ",
             "'mirnaObj'. Please, use 'performGeneDE()' before using ",
             "this function. See ?performGeneDE",
             call. = FALSE
        )
    }
    if (max(dim(mirnaTargets(mirnaObj))) == 0) {
        stop("No targets are present within 'mirnaObj'! ",
             "Before performing the integration analysis miRNA target ",
             "genes must be retrieved with the 'getTargets()' function. ",
             "See '?getTargets' for the details.",
             call. = FALSE
        )
    }
    if (!is.character(test) |
        length(test) != 1 |
        !test %in% c("auto", "correlation", "association", "fry")) {
        stop("'test' must be one of:\n",
             "\t- 'auto', (default) to automatically ",
             "choose the appropriate test;\n",
             "\t- 'correlation', to perform a correlation analysis;\n",
             "\t- 'association', to apply a one-sided association test;\n",
             "\t- 'fry', to apply rotation gene-set testing.",
             call. = FALSE
        )
    }
    if (geneDE(mirnaObj, param = TRUE)$method == "Manually added" &
        test == "fry") {
        stop("Integration through 'fry' method is not available ",
             "for user-supplied differential expression results...",
             call. = FALSE
        )
    }
    if (test == "correlation" & pairedSamples(mirnaObj) == FALSE) {
        warning("You can't perform a 'correlation' analysis with ",
                "unpaired samples! Setting 'test' to 'association' for ",
                "performing a one-sided association test instead...",
                call. = FALSE
        )
        test <- "association"
    }
    if (!is.numeric(pCutoff) |
        length(pCutoff) != 1 |
        pCutoff > 1 |
        pCutoff < 0) {
        stop("'pCutoff' must be a number between 0 and 1! (default is 0.05)",
             call. = FALSE
        )
    }
    if (!is.character(pAdjustment) |
        length(pAdjustment) != 1 |
        !pAdjustment %in% c(
            "none", "fdr", "bonferroni", "BY", "hochberg",
            "holm", "hommel", "BH"
        )) {
        stop("'pAdjustment' must be  one of: 'none', 'fdr' (default), ",
             "'BH' (same as 'fdr'), 'bonferroni', 'BY', 'hochberg', ",
             "'holm', 'hommel'",
             call. = FALSE
        )
    }
    if (!is.character(corMethod) |
        length(corMethod) != 1 |
        !corMethod %in% c("pearson", "kendall", "spearman")) {
        stop("'corMethod' must be one of: 'pearson', 'kendall' and ",
             "'spearman' (default).",
             call. = FALSE
        )
    }
    if (!is.numeric(corCutoff) |
        length(corCutoff) != 1 |
        corCutoff > 1 |
        corCutoff < 0) {
        stop("'corCutoff' must be a number between 0 and 1! (default is 0.5)",
             call. = FALSE
        )
    }
    if (!is.character(associationMethod) |
        length(associationMethod) != 1 |
        !associationMethod %in% c("boschloo", "fisher-midp", "fisher")) {
        stop("'associationMethod' must be one of: 'boschloo' (default), ",
             "'fisher-midp', and 'fisher'.",
             call. = FALSE
        )
    }
    if (!is.numeric(nuisanceParam) |
        length(nuisanceParam) != 1 |
        nuisanceParam < 0) {
        stop("'nuisanceParam' must be a non negative number! (default is 100)",
             call. = FALSE
        )
    }
    
    ## use the appropriate test
    if (test == "auto" & pairedSamples(mirnaObj) == TRUE) {
        intMethod <- "correlation"
        message(
            "Since data derive from paired samples, a correlation test ",
            "will be used."
        )
    } else if (test == "auto" & pairedSamples(mirnaObj) == FALSE) {
        intMethod <- "association"
        message(
            "Since data derive from different individuals, a ",
            "one-sided association test will be used."
        )
    } else {
        intMethod <- test
        message(
            "As specified by the user, a ", test,
            " will be used."
        )
    }
    
    ## call the right function
    if (intMethod == "correlation") {
        mirnaObj <- correlateMirnaTargets(
            mirnaObj,
            corMethod,
            corCutoff,
            pCutoff,
            pAdjustment
        )
    } else if (intMethod == "association") {
        mirnaObj <- associateMirnaTargets(
            mirnaObj,
            pCutoff,
            pAdjustment,
            associationMethod,
            nuisanceParam
        )
    } else if (intMethod == "fry") {
        mirnaObj <- fryMirnaTargets(
            mirnaObj,
            pCutoff,
            pAdjustment
        )
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
    sMap <- sampleMap(mirnaObj)
    mirnaSamples <- sMap$primary[sMap$assay == "microRNA"]
    geneSamples <- sMap$primary[sMap$assay == "genes"]
    
    if (!identical(mirnaSamples, geneSamples)) {
        ## determine common and uncommon samples
        common <- intersect(mirnaSamples, geneSamples)
        unpaired <- setdiff(c(mirnaSamples, geneSamples), common)
        
        ## stop if common samples are less than 3
        if (length(common) < 3) {
            stop("To perform the correlation at least 3 paired samples ",
                 "are needed.",
                 call. = FALSE
            )
        }
        
        ## remove samples without measurments of both miRNAs and genes
        if (length(unpaired) > 0) {
            mirnaExpr <- mirnaExpr[, sMap$colname[sMap$assay == "microRNA" &
                                                      sMap$primary %in% common]]
            geneExpr <- geneExpr[, sMap$colname[sMap$assay == "genes" &
                                                    sMap$primary %in% common]]
            
            warning("Some samples don't have expression values for both ",
                    "miRNAs and genes, thus they have been excluded from the ",
                    "correlation analysis. Removed samples are: ",
                    paste(unpaired, collapse = ", "),
                    call. = FALSE
            )
        }
        
        ## order the columns of expression matrices in the same way
        mirnaMap <- sMap[sMap$assay == "microRNA" &
                             sMap$primary %in% common, ]
        geneMap <- sMap[sMap$assay == "genes" &
                            sMap$primary %in% common, ]
        mirnaOrder <- mirnaMap$primary[order(match(
            mirnaMap$colname,
            colnames(mirnaExpr)
        ))]
        geneExpr <- geneExpr[, geneMap$colname[order(match(
            geneMap$primary,
            mirnaOrder
        ))]]
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
    message("Performing ", usedCoef, "'s correlation analysis...")
    correlation <- mapply(function(mirna, gene) {
        ## extract the expression values of the microRNA and the target
        mirnaInt <- as.numeric(mirnaExpr[mirna, ])
        geneInt <- as.numeric(targetExpr[gene, ])
        
        ## perform the correlation analysis
        corPair <- cor.test(mirnaInt,
                            geneInt,
                            method = corMethod,
                            alternative = "less",
                            exact = FALSE
        )
        
        ## report the results of the correlation analysis
        fold <- ifelse(dem$logFC[dem$ID == mirna] > 0,
                       "upregulated",
                       "downregulated"
        )
        
        pair <- c(
            mirna,
            gene,
            fold,
            corPair$estimate,
            corPair$p.value
        )
        pair
    }, targetsTable$MicroRNA, targetsTable$Gene.Symbol)
    
    ## convert correlation output to a data.frame object
    corRes <- as.data.frame(t(correlation))
    colnames(corRes) <- c(
        "microRNA", "Target", "microRNA.Direction",
        "Corr.Coefficient", "Corr.P.Value"
    )
    corRes$Corr.Coefficient <- as.numeric(corRes$Corr.Coefficient)
    corRes$Corr.P.Value <- as.numeric(corRes$Corr.P.Value)
    
    ## correct correlation p values for multiple testing
    pAdj <- p.adjust(corRes$Corr.P.Value, method = pAdjustment)
    corRes$Corr.Adjusted.P.Val <- pAdj
    
    ## select statistically significant associations
    corRes <- corRes[corRes$Corr.Adjusted.P.Val <= pCutoff &
                         abs(corRes$Corr.Coefficient) >= corCutoff, ]
    
    ## report the results of the correlation analysis
    if (nrow(corRes) >= 1) {
        message(
            "A statistically significant correlation between ",
            nrow(corRes), " miRNA-target pairs was found!"
        )
    } else {
        message(
            "No statistically significant correlation between ",
            "DE-miRNAs and targets was found."
        )
    }
    
    ## return the object with the results of the correlation analysis
    resList <- list(
        data = corRes,
        method = paste(usedCoef, "'s correlation analysis", sep = ""),
        pCutoff = pCutoff,
        pAdjustment = pAdjustment
    )
    integration(mirnaObj) <- resList
    return(mirnaObj)
}





## one-sided association tests for integrative analysis
associateMirnaTargets <- function(mirnaObj,
                                  pCutoff,
                                  pAdjustment,
                                  associationMethod,
                                  nuisanceParam) {
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
    if (associationMethod == "boschloo") {
        meth <- "One-sided Boschloo's exact test"
    } else if (associationMethod == "fisher-midp") {
        meth <- paste(
            "One-sided Fisher's exact test",
            "with Lancaster's mid-p correction"
        )
    } else if (associationMethod == "fisher") {
        meth <- "One-sided Fisher's exact test"
    }
    
    ## compute the inverse association between each miRNA and its targets
    message("Performing ", meth, "...")
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
        
        contingencyTable <- matrix(c(
            nonDegsTarg,
            nonDegsNonTarg,
            degsTarg,
            degsNonTarg
        ), 2)
        
        ## perform the appropriate statistical test
        pval <- association.helper(
            contingencyTable,
            associationMethod,
            nuisanceParam
        )
        
        ## report results to a list object
        intRes <- c(
            x,
            fold,
            geneFold,
            degsTarg,
            degsTarg + nonDegsTarg,
            pval,
            paste(intersect(degs, targ), collapse = "/")
        )
        intRes
    })
    
    ## convert list object into a data.frame
    association <- as.data.frame(do.call(rbind, association))
    colnames(association) <- c(
        "microRNA",
        "mirna.direction",
        "gene.direction",
        "DE",
        "targets",
        "P.Val",
        "DE.targets"
    )
    
    ## convert numeric columns
    association$P.Val <- as.numeric(association$P.Val)
    association$DE <- as.numeric(association$DE)
    association$targets <- as.numeric(association$targets)
    
    ## correct p-values
    pAdj <- p.adjust(association$P.Val, method = pAdjustment)
    association$adj.P.Val <- pAdj
    association <- association[, c(1, 2, 3, 4, 5, 6, 8, 7)]
    association <- association[order(association$adj.P.Val), ]
    
    ## print association results
    association <- association[association$adj.P.Val <= pCutoff, ]
    if (nrow(association) >= 1) {
        message(
            "A statistically significant association between ",
            nrow(association), " DE-miRNAs and ",
            length(unique(unlist(strsplit(association$DE.targets, "/")))),
            " genes was found!"
        )
    } else {
        message(
            "No statistically significant associations between ",
            "DE-miRNAs and DEGs were found."
        )
    }
    
    ## return the object containing association results
    resList <- list(
        data = association,
        method = meth,
        pCutoff = pCutoff,
        pAdjustment = pAdjustment
    )
    integration(mirnaObj) <- resList
    return(mirnaObj)
}





## helper function for choosing the association test
association.helper <- function(
        contingencyTable,
        associationMethod,
        nuisanceParam) {
    ## perform the appropriate statistical test
    if (associationMethod == "boschloo") {
        pval <- boshloo.test(contingencyTable, nuisanceParam)
    } else if (associationMethod == "fisher-midp") {
        pval <- fisher.midp(contingencyTable, midpAdjustment = TRUE)
    } else if (associationMethod == "fisher") {
        pval <- fisher.midp(contingencyTable, midpAdjustment = FALSE)
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
        midp <- phyper(mat, m, n, k) - 0.5 * dhyper(mat, m, n, k)
    } else {
        midp <- phyper(mat, m, n, k)
    }
    return(midp)
}





## optimization helper function for p-value refinement
optimization.helper <- function(p, Ns, moreExtremeMat) {
    sum(matrix(dbinom(seq(0, Ns[1]), Ns[1], p), ncol = 1) *
            (moreExtremeMat %*% dbinom(seq(0, Ns[2]), Ns[2], p)))
}





## calculate p-values fro Boshloo's exact test
boshloo.test <- function(data, npNumbers) {
    ## this implementation of the Boschloo's test is
    ## inspired by the Exact R package (Peter Calhoun)
    
    ## create a vector of values for the nusiance parameter
    int <- seq(0.00001, 0.99999, length = npNumbers)
    
    ## transpose contingency table to have fixed row sums
    data <- t(data)
    
    ## extract values from contingency table
    x <- data[1, 1]
    y <- data[1, 2]
    Ns <- .colSums(data, 2, 2)
    N <- sum(Ns)
    
    ## compute p-value resulting from Fisher's exact test
    fp <- fisher.midp(data, midpAdjustment = FALSE)
    
    ## define a matrix for calculating the number of more extreme tables
    moreExtremeMat <- matrix(NA, Ns[1] + 1, Ns[2] + 1,
                             dimnames = list(seq(0, Ns[1]), seq(0, Ns[2]))
    )
    for (j in seq(data[1, 2] + 1, Ns[2] + 1)) {
        moreExtremeMat[seq(0, data[1, 1] + 1), j] <- 1
    }
    
    ## compute the number of more extreme tables for each row of the matrix
    for (i in seq(0, Ns[1])) {
        ## find first column with NA
        startJ <- which(is.na(moreExtremeMat[i + 1, ]))[1] - 1
        if (!is.na(startJ)) {
            ## compute the number of more extreme tables for each column
            for (j in startJ:Ns[2]) {
                ## define a new contingency table for each table
                newDat <- matrix(c(i, Ns[1] - i, j, Ns[2] - j), 2, 2)
                
                ## calculate Fisher's exact test p-value for each table
                newTX <- fisher.midp(newDat, midpAdjustment = FALSE)
                newTX <- signif(newTX, 12) ## correct for rounding errors
                
                ## compare Fisher's p-value to the observed p-value
                if (newTX <= fp) {
                    ## set the remaining columns in the row as more extreme
                    moreExtremeMat[i + 1, seq(j + 1, Ns[2] + 1)] <- 1
                    break
                } else {
                    ## set the remaining rows in the column as less extreme
                    moreExtremeMat[seq(i + 1, Ns[1] + 1), j + 1] <- 0
                }
            }
        }
    }
    
    ## create probability vector
    index <- 1
    prob <- rep(NA, length(int))
    
    ## select the binomials needed
    Tbls <- which(moreExtremeMat == 1, arr.ind = TRUE) - 1
    maxX1 <- max(Tbls[, 1])
    minX2 <- min(Tbls[, 2])
    
    ## calculate the probability for each value of the nuisance parameter
    for (probVal in int) {
        binomProb1 <- dbinom(seq(0, maxX1), Ns[1], probVal)
        binomProb2 <- dbinom(seq(minX2, Ns[2]), Ns[2], probVal)
        matVal <- moreExtremeMat[seq(maxX1 + 1), seq(minX2 + 1, Ns[2] + 1),
                                 drop = FALSE]
        prob[index] <- suppressWarnings(
            sum(binomProb1 * (matVal %*% binomProb2)))
        index <- index + 1
    }
    prob <- signif(prob, 12)
    
    ## identify the maximum p-value
    np <- int[which(prob == max(prob, na.rm = TRUE))]
    pvalue <- max(prob, na.rm = TRUE)
    
    ## refine p-values through the optimize function
    refPvalue <- rep(0, length(np))
    refNp <- refPvalue
    for (i in seq_along(np)) {
        ref <- suppressWarnings(
            optimize(
                f = optimization.helper,
                interval = c(
                    max(int[1], np[i] - 1 / npNumbers),
                    min(int[npNumbers], np[i] + 1 / npNumbers)
                ),
                Ns = Ns,
                moreExtremeMat = moreExtremeMat,
                maximum = TRUE
            )
        )
        refPvalue[i] <- ref$objective
        refNp[i] <- ref$maximum
    }
    refPvalue <- min(c(1, signif(refPvalue, 12)))
    if (!all(is.na(refPvalue)) && max(refPvalue, na.rm = TRUE) > pvalue) {
        np <- refNp[refPvalue == max(refPvalue)]
        pvalue <- max(refPvalue)
    }
    
    ## return p-value
    return(pvalue)
}





## rotation gene-set test for integrative analysis
fryMirnaTargets <- function(mirnaObj,
                            pCutoff,
                            pAdjustment) {
    ## extract gene differential expression results
    de <- geneDE(mirnaObj, param = TRUE)
    
    ## access sample metadata
    meta <- colData(mirnaObj)
    meta <- meta[!is.na(meta$geneCol), ]
    
    ## reorder metadata based on expression matrix
    meta <- meta[order(match(meta$geneCol, colnames(mirnaObj[["genes"]]))), ]
    
    ## build miRNA-target sets
    tg <- mirnaTargets(mirnaObj)
    tgList <- split(tg$Gene.Symbol, tg$MicroRNA)
    
    ## determine the appropriate expression matrix and the experimental design
    if (de$method == "limma") {
        expr <- mirnaObj[["genes"]]
        des <- model.matrix(de$design, data = meta)
    } else if (de$method == "edgeR" |
               de$method == "limma-voom") {
        expr <- geneDE(mirnaObj, returnObject = TRUE)
        des <- expr$design
    } else if (de$method == "DESeq2") {
        message("Applying 'limma-voom' pipeline before using 'fry' method...")
        counts <- metadata(mirnaObj)[["oldCounts"]][["genes"]]
        features <- edgeR::DGEList(
            counts = counts,
            group = meta[, de$group],
            samples = meta
        )
        keep <- edgeR::filterByExpr(features)
        features <- features[keep, , keep.lib.sizes = FALSE]
        features <- edgeR::calcNormFactors(features)
        des <- model.matrix(de$design, data = meta)
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
        con <- limma::makeContrasts(
            contrasts = contrast,
            levels = des
        )
    } else {
        ## set the coefficient name for the appropriate comparison
        con <- contrast[1]
        con <- paste(de$group, con, sep = "")
    }
    
    ## perform integration through rotation gene set enrichment analysis
    message("Performing miRNA-gene integrative analysis using 'fry' method...")
    if (de$method == "edgeR") {
        ## perform miRNA-target integration through 'fry'
        rs <- edgeR::fry.DGEList(
            y = expr,
            index = tgList,
            design = des,
            contrast = con,
            adjust.method = pAdjustment
        )
    } else {
        ## perform miRNA-target integration through 'fry'
        rs <- limma::fry(
            y = expr,
            index = tgList,
            design = des,
            contrast = con,
            adjust.method = pAdjustment
        )
    }
    
    ## retain effects in the right direction
    dem <- mirnaDE(mirnaObj)
    upDem <- dem$ID[dem$logFC > 0]
    downDem <- dem$ID[dem$logFC < 0]
    rs <- rs[(rownames(rs) %in% upDem & rs$Direction == "Down") |
                 (rownames(rs) %in% downDem & rs$Direction == "Up"), ]
    
    ## maintain interactions under the specified cutoff
    res <- rs[rs$FDR <= pCutoff, ]
    
    ## print integration results
    if (nrow(res) == 0) {
        res <- data.frame()
        message(
            "No statistically significant associations between ",
            "DE-miRNAs and genes were found."
        )
    } else {
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
        colnames(res) <- c(
            "microRNA", "mirna.direction", "gene.direction",
            "DE", "targets", "P.Val", "adj.P.Val", "DE.targets"
        )
        
        ## report results
        message(
            "A statistically significant association between ",
            nrow(res), " DE-miRNAs and ",
            length(unique(unlist(strsplit(res$DE.targets, "/")))),
            " genes was found!"
        )
    }
    
    ## return the object containing association results
    resList <- list(
        data = res,
        method = "Rotation Gene-Set Test (FRY)",
        pCutoff = pCutoff,
        pAdjustment = pAdjustment
    )
    integration(mirnaObj) <- resList
    return(mirnaObj)
}
