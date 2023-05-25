#' List the available miEAA 2.0 categories for an organism
#'
#' This function retrieves a list of miEAA 2.0 categories that are available
#' for a specified organism. The name of these categories can be later used
#' with [enrichMirnas()] and [gseaMirnas()] functions to perform miRNA
#' enrichment analysis with the desired category.
#'
#' @param organism The name of the organism under consideration. It must be one
#' of: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`,
#' `Arabidopsis thaliana`, `Bos taurus`, `Caenorhabditis elegans`,
#' `Drosophila melanogaster`, `Danio rerio`, `Gallus gallus`, `Sus scrofa`
#' @param mirnaType The type of miRNA specie: use `mature` for the enrichment
#' of mature microRNAs; use `precursor`, for the enrichment analysis of miRNA
#' precursors
#'
#' @returns
#' A named vector containing the IDs of miEAA 2.0 categories. The user can
#' then decide the categories to use with [enrichMirnas()] and [gseaMirnas()]
#' functions.
#'
#' @examples
#' # list all the available categories for Mus musculus
#' available_categories <- listCategories(organism = "Mus musculus")
#'
#' @note
#' This function uses the package `rbioapi` to interface miEAA 2.0 database
#' for retrieving the available categories for different organisms.
#'
#' @references
#' Moosa Rezwani, Ali Akbar Pourfathollah, Farshid Noorbakhsh, rbioapi:
#' user-friendly R interface to biologic web services’ API, Bioinformatics,
#' Volume 38, Issue 10, 15 May 2022, Pages 2952–2953,
#' \url{https://doi.org/10.1093/bioinformatics/btac172}
#'
#' Fabian Kern, Tobias Fehlmann, Jeffrey Solomon, Louisa Schwed, Nadja Grammes,
#' Christina Backes, Kendall Van Keuren-Jensen, David Wesley Craig, Eckart
#' Meese, Andreas Keller, miEAA 2.0: integrating multi-species microRNA
#' enrichment analysis and workflow management systems, Nucleic Acids Research,
#' Volume 48, Issue W1, 02 July 2020, Pages W521–W528,
#' \url{https://doi.org/10.1093/nar/gkaa309}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
listCategories <- function(organism = "Homo sapiens",
                           mirnaType = "mature") {

  ## check inputs
  if (!organism %in% convertOrganism("miEAA", "all")) {
    stop(paste("'organism' must be one of:",
               paste(convertOrganism("miEAA", "all"), collapse = ", ")),
         call. = FALSE)
  }
  if (!is.character(mirnaType) |
      length(mirnaType) != 1 |
      !mirnaType %in% c("mature", "precursor")) {
    stop("'mirnaType' must be one of: 'mature', 'precursor'.", call. = FALSE)
  }

  ## retrieve miEAA categories
  cats <- rbioapi::rba_mieaa_cats(mirna_type = mirnaType,
                                  species = organism)

  ## return categories
  return(cats)

}





#' Perform an over-representation analysis (ORA) for DE-miRNAs
#'
#' This function uses miEAA 2.0 (Kern et al., 2020) to perform an
#' over-representation analysis of differentially expressed miRNAs in order to
#' investigate the biological consequences of miRNA dysregulations. This can be
#' pivotal for understanding the biological pathways affected by
#' miRNA alterations.
#'
#' @details
#' To understand the functions of DE-miRNAs, one of the most common approaches
#' is to retrieve DE-miRNA targets and perform the enrichment of these genes.
#' However, Bleazard et al. discovered that this approach may lead to biased
#' results, thus limiting the reliability of the enrichment analysis.
#' To overcome this limit, several other resources were developed, including
#' miEAA 2.0, which is a tool based on GeneTrail that facilitates the
#' interpretation of microRNA dysregulations.
#' This tool allows to perform enrichment analyses with multiple categories for
#' assessing different biological features.
#' These categories, which may not be available for some species,
#' include databases like miRWalk, mirPathDB and MNDR, and allow to investigate
#' different features such as biological processes, molecular functions,
#' biological pathways, diseases, target genes, organs, tissues and more. To
#' retrieve the complete list of available miEAA 2.0 categories, use the
#' [listCategories()] function.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param organism The name of the organism under consideration. It must be one
#' of: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`,
#' `Arabidopsis thaliana`, `Bos taurus`, `Caenorhabditis elegans`,
#' `Drosophila melanogaster`, `Danio rerio`, `Gallus gallus`, `Sus scrofa`
#' @param category The name of the miEAA 2.0 category to use for the enrichment
#' analysis. The list of available categories can be obtained through the
#' [listCategories()] function. Moreover, the user can also perform the
#' enrichment analysis on all miEAA 2.0 categories by specifying `all` (default)
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param independentAdj Logical, whether to adjust the p-values of each
#' category independently or to adjust p-values collectively over all
#' categories. Default is `TRUE`
#' @param minHits The minimum number of miRNAs that must be present in a
#' subcategory in order to include that subcategory in the results
#' (default is `2`)
#' @param mirnaType The type of miRNA specie: use `mature` for the enrichment
#' of mature microRNAs; use `precursor`, for the enrichment analysis of miRNA
#' precursors
#' @param sortBy The column used to sort the results table. It must be one of:
#' `category` (default), `subcategory`, `enrichment`, `p_value`, `p_adjusted`,
#' `q_value`, `observed`
#' @param convertVer Logical, whether miRNA IDs stored in `mirnaObj` need to
#' be converted to miRBase version 22. Default is `FALSE`
#' @param mirBaseVer The miRBase version of miRNA IDs stored in `mirnaObj`.
#' This parameter is only needed if miRNA IDs need to be converted to
#' version 22
#' @param ... Other arguments that can be passed to
#' [rbioapi::rba_mieaa_enrich()]
#'
#' @returns
#' A `list` object containing miEAA 2.0 enrichment results. This `list` object
#' has two [`MirnaEnrichment`][MirnaEnrichment-class] objects, one for the
#' enrichment results of up-regulated miRNAs and one for enrichment results of
#' down-regulated miRNAs. To access one or the other, the user can make use of
#' `res[["upregulated"]]` and `res[["downregulated"]]`. Further, it is possible
#' to access the results of these [`MirnaEnrichment`][MirnaEnrichment-class]
#' objects through the [enrichmentResults()] function. For additional
#' instructions see [`MirnaEnrichment`][MirnaEnrichment-class] class.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform an over-represantation analysis with miEAA category GO Annotations
#' enr <- enrichMirnas(obj, organism = "Homo sapiens",
#' category = "GO_Annotations_mature")
#'
#' # divide the enrichment of up- and down-regulated miRNAs
#' enr_up <- enr[["upregulated"]]
#' enr_down <- enr[["downregulated"]]
#'
#' # extract results
#' res_up <- enrichmentResults(enr_up)
#'
#' # plot results
#' mirnaDotplot(enr_up)
#'
#' @note
#' To query miEAA 2.0, this function uses the package `rbioapi`, which
#' provides convenient methods for accessing miEAA 2.0 to perform
#' enrichment analyses.
#'
#' @references
#' Fabian Kern, Tobias Fehlmann, Jeffrey Solomon, Louisa Schwed, Nadja Grammes,
#' Christina Backes, Kendall Van Keuren-Jensen, David Wesley Craig, Eckart
#' Meese, Andreas Keller, miEAA 2.0: integrating multi-species microRNA
#' enrichment analysis and workflow management systems, Nucleic Acids Research,
#' Volume 48, Issue W1, 02 July 2020, Pages W521–W528,
#' \url{https://doi.org/10.1093/nar/gkaa309}
#'
#' Thomas Bleazard, Janine A Lamb, Sam Griffiths-Jones, Bias in microRNA
#' functional enrichment analysis, Bioinformatics, Volume 31, Issue 10, 15 May
#' 2015, Pages 1592–1598, \url{https://doi.org/10.1093/bioinformatics/btv023}
#'
#' Moosa Rezwani, Ali Akbar Pourfathollah, Farshid Noorbakhsh, rbioapi:
#' user-friendly R interface to biologic web services’ API, Bioinformatics,
#' Volume 38, Issue 10, 15 May 2022, Pages 2952–2953,
#' \url{https://doi.org/10.1093/bioinformatics/btac172}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichMirnas <- function(mirnaObj,
                         organism = "Homo sapiens",
                         category = "all",
                         pCutoff = 0.05,
                         pAdjustment = "fdr",
                         independentAdj = TRUE,
                         minHits = 2,
                         mirnaType = "mature",
                         sortBy = "category",
                         convertVer = FALSE,
                         mirBaseVer = 22,
                         ...) {

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
  if (!organism %in% convertOrganism("miEAA", "all")) {
    stop(paste("'organism' must be one of:",
               paste(convertOrganism("miEAA", "all"), collapse = ", ")),
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
                          "holm", "hommel")) {
    stop(paste("'pAdjustment' must be  one of: 'none', 'fdr' (default),",
               "'bonferroni', 'BY', 'hochberg', 'holm', 'hommel'"),
         call. = FALSE)
  }
  if (!is.logical(independentAdj) |
      length(independentAdj) != 1) {
    stop("'independentAdj' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(minHits) |
      length(minHits) != 1) {
    stop("'minHits' must be a numeric object of lenght one (default is 2)",
         call. = FALSE)
  }
  if (!is.character(mirnaType) |
      length(mirnaType) != 1 |
      !mirnaType %in% c("mature", "precursor")) {
    stop("'mirnaType' must be one of: 'mature', 'precursor'.", call. = FALSE)
  }
  if (!is.character(sortBy) |
      length(sortBy) != 1 |
      !sortBy %in% c("category", "subcategory", "enrichment", "p_value",
                     "p_adjusted", "q_value", "observed")) {
    stop(paste("'sortBy' must be one of: 'category', 'subcategory',",
               "'enrichment', 'p_value', 'p_adjusted',",
               "'q_value', 'observed'."), call. = FALSE)
  }
  if (!is.logical(convertVer) |
      length(convertVer) != 1) {
    stop("'convertVer' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(mirBaseVer) |
      length(mirBaseVer) != 1) {
    stop(paste("'mirBaseVer' must be a numeric object of lenght one",
               "specifying the miRBase version of miRNA IDs in 'mirnaObj'."),
         call. = FALSE)
  }

  validCategories <- suppressMessages(
    rbioapi::rba_mieaa_cats(mirna_type = mirnaType,
                            species = organism)
  )
  if (!category %in% validCategories & category != "all") {
    stop(paste("'category' must be one of:\n", "\n- all",
               paste(paste("-", validCategories), collapse = "\n"), sep = ""),
         call. = FALSE)
  }

  if (category == "all") {
    cat <- NULL
  } else {
    cat <- category
  }

  ## extract DE-miRNAs
  mirna <- mirnaDE(mirnaObj)
  upMirna <- mirna[mirna$logFC > 0, "ID"]
  downMirna <- mirna[mirna$logFC < 0, "ID"]

  ## convert miRNAs to miRBase v22
  if (convertVer == TRUE) {
    upMirna <- quiet(suppressMessages(
      rbioapi::rba_mieaa_convert_version(upMirna,
                                         mirna_type = mirnaType,
                                         input_version = mirBaseVer,
                                         output_version = 22,
                                         simple_output = TRUE,
                                         ...)
    ))
    downMirna <- quiet(suppressMessages(
      rbioapi::rba_mieaa_convert_version(downMirna,
                                         mirna_type = mirnaType,
                                         input_version = mirBaseVer,
                                         output_version = 22,
                                         simple_output = TRUE,
                                         ...)
    ))
  }

  ## set background genes
  universe <- rownames(mirnaObj)[["microRNA"]]

  ## retrieve miRNA enrichment
  message("Querying miEAA 2.0 for enriching up-regulated microRNAs\n")
  enrUp <- quiet(suppressMessages(
    rbioapi::rba_mieaa_enrich(test_set = upMirna,
                              ref_set = universe,
                              mirna_type = mirnaType,
                              test_type = "ORA",
                              species = organism,
                              categories = cat,
                              p_adj_method = pAdjustment,
                              independent_p_adj = independentAdj,
                              sig_level = pCutoff,
                              min_hits = minHits,
                              sort_by = sortBy,
                              ...)
  ))
  message("Querying miEAA 2.0 for enriching down-regulated microRNAs\n")
  enrDown <- quiet(suppressMessages(
    rbioapi::rba_mieaa_enrich(test_set = downMirna,
                              ref_set = universe,
                              mirna_type = mirnaType,
                              test_type = "ORA",
                              species = organism,
                              categories = cat,
                              p_adj_method = pAdjustment,
                              independent_p_adj = independentAdj,
                              sig_level = pCutoff,
                              min_hits = minHits,
                              sort_by = sortBy,
                              ...)
  ))

  ## convert wrong character columns to numeric
  numCols <- c("P-value", "P-adjusted", "Q-value", "Expected", "Observed")
  if (nrow(enrUp) > 0) {
    enrUp[, numCols] <- sapply(enrUp[, numCols], as.numeric)
  }
  if (nrow(enrDown) > 0) {
    enrDown[, numCols] <- sapply(enrDown[, numCols], as.numeric)
  }

  ## rename columns with consistent names
  colnames(enrUp) <- gsub("-", ".", colnames(enrUp))
  colnames(enrDown) <- gsub("-", ".", colnames(enrDown))

  ## creation of output objects
  upRes <- new("MirnaEnrichment",
               result = enrUp,
               pvalueCutoff = pCutoff,
               pAdjustMethod = pAdjustment,
               organism = organism,
               ontology = category,
               gene = upMirna,
               universe = universe,
               keytype = "miRBase v22")
  downRes <- new("MirnaEnrichment",
                 result = enrDown,
                 pvalueCutoff = pCutoff,
                 pAdjustMethod = pAdjustment,
                 organism = organism,
                 ontology = category,
                 gene = downMirna,
                 universe = universe,
                 keytype = "miRBase v22")

  ## print a message with the results of the enrichment
  message(paste("\nThe functional enrichment analysis of miRNAs reported",
                nrow(enrDown),
                "significantly enriched terms for downregulated miRNAs and",
                nrow(enrUp),
                "for upregulated miRNAs.\n"))

  ## return a list object with enrichments of up- and downregulated miRNAs
  res <- list("upregulated" = upRes, "downregulated" = downRes)
  return(res)

}





#' Perform a gene set enrichment analysis (GSEA) with miRNA expression values
#'
#' This function uses miEAA 2.0 (Kern et al., 2020) to perform a gene set
#' enrichment analysis (GSEA) with miRNA expression data in order to
#' investigate the biological consequences of miRNA dysregulations. This can be
#' pivotal for understanding the biological pathways affected by
#' miRNA alterations.
#'
#' @details
#' Recently, gene set enrichment analysis (GSEA) has emerged to be a valuable
#' tool for assessing the biological effects of gene expression dysregulations.
#' However, few resources implemented the possibility of running this algorithm
#' with miRNA expression data. One tool that provides this feature is miEAA 2.0,
#' which adopts an un-weighted variant of the original algorithm, with constant
#' running sum changes, resulting in a Kolmogorow-Smirnow test.
#' miEAA 2.0, allows to perform enrichment analyses with multiple categories
#' for assessing different biological features. These categories, which may not
#' be available for some species, include databases like miRWalk, mirPathDB and
#' MNDR, and allow to investigate different features such as biological
#' processes, molecular functions, biological pathways, diseases, target genes,
#' organs, tissues and more. To retrieve the complete list of available
#' miEAA 2.0 categories, use the [listCategories()] function.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param organism The name of the organism under consideration. It must be one
#' of: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`,
#' `Arabidopsis thaliana`, `Bos taurus`, `Caenorhabditis elegans`,
#' `Drosophila melanogaster`, `Danio rerio`, `Gallus gallus`, `Sus scrofa`
#' @param category The name of the miEAA 2.0 category to use for the enrichment
#' analysis. The list of available categories can be obtained through the
#' [listCategories()] function. Moreover, the user can also perform the
#' enrichment analysis on all miEAA 2.0 categories by specifying `all` (default)
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param independentAdj Logical, whether to adjust the p-values of each
#' category independently or to adjust p-values collectively over all
#' categories. Default is `TRUE`
#' @param minHits The minimum number of miRNAs that must be present in a
#' subcategory in order to include that subcategory in the results
#' (default is `2`)
#' @param mirnaType The type of miRNA specie: use `mature` for the enrichment
#' of mature microRNAs; use `precursor`, for the enrichment analysis of miRNA
#' precursors
#' @param sortBy The column used to sort the results table. It must be one of:
#' `category` (default), `subcategory`, `enrichment`, `p_value`, `p_adjusted`,
#' `q_value`, `observed`
#' @param convertVer Logical, whether miRNA IDs stored in `mirnaObj` need to
#' be converted to miRBase version 22. Default is `FALSE`
#' @param mirBaseVer The miRBase version of miRNA IDs stored in `mirnaObj`.
#' This parameter is only needed if miRNA IDs need to be converted to
#' version 22
#' @param ... Other arguments that can be passed to
#' [rbioapi::rba_mieaa_enrich()]
#'
#' @returns
#' A [`MirnaGsea`][MirnaGsea-class] object containing the results of gene set
#' enrichment analysis. To access the results table with the enriched terms,
#' the user can make use of the [enrichmentResults()] function.
#' For additional instructions see [`MirnaGsea`][MirnaGsea-class] class.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform a GSEA with the miEAA category GO Annotations
#' gse <- gseaMirnas(obj, organism = "Homo sapiens",
#' category = "GO_Annotations_mature")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' mirnaDotplot(gse)
#' mirnaRidgeplot(gse)
#'
#' @note
#' To query miEAA 2.0, this function uses the package `rbioapi`, which
#' provides convenient methods for accessing miEAA 2.0 to perform
#' enrichment analyses.
#'
#' @references
#' Fabian Kern, Tobias Fehlmann, Jeffrey Solomon, Louisa Schwed, Nadja Grammes,
#' Christina Backes, Kendall Van Keuren-Jensen, David Wesley Craig, Eckart
#' Meese, Andreas Keller, miEAA 2.0: integrating multi-species microRNA
#' enrichment analysis and workflow management systems, Nucleic Acids Research,
#' Volume 48, Issue W1, 02 July 2020, Pages W521–W528,
#' \url{https://doi.org/10.1093/nar/gkaa309}
#'
#' Thomas Bleazard, Janine A Lamb, Sam Griffiths-Jones, Bias in microRNA
#' functional enrichment analysis, Bioinformatics, Volume 31, Issue 10, 15 May
#' 2015, Pages 1592–1598, \url{https://doi.org/10.1093/bioinformatics/btv023}
#'
#' Moosa Rezwani, Ali Akbar Pourfathollah, Farshid Noorbakhsh, rbioapi:
#' user-friendly R interface to biologic web services’ API, Bioinformatics,
#' Volume 38, Issue 10, 15 May 2022, Pages 2952–2953,
#' \url{https://doi.org/10.1093/bioinformatics/btac172}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
gseaMirnas <- function(mirnaObj,
                       organism = "Homo sapiens",
                       category = "all",
                       pCutoff = 0.05,
                       pAdjustment = "fdr",
                       independentAdj = TRUE,
                       minHits = 2,
                       mirnaType = "mature",
                       sortBy = "category",
                       convertVer = FALSE,
                       mirBaseVer = 22,
                       ...) {

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
  if (!organism %in% convertOrganism("miEAA", "all")) {
    stop(paste("'organism' must be one of:",
               paste(convertOrganism("miEAA", "all"), collapse = ", ")),
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
                          "holm", "hommel")) {
    stop(paste("'pAdjustment' must be  one of: 'none', 'fdr' (default),",
               "'bonferroni', 'BY', 'hochberg', 'holm', 'hommel'"),
         call. = FALSE)
  }
  if (!is.logical(independentAdj) |
      length(independentAdj) != 1) {
    stop("'independentAdj' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(minHits) |
      length(minHits) != 1) {
    stop("'minHits' must be a numeric object of lenght one (default is 2)",
         call. = FALSE)
  }
  if (!is.character(mirnaType) |
      length(mirnaType) != 1 |
      !mirnaType %in% c("mature", "precursor")) {
    stop("'mirnaType' must be one of: 'mature', 'precursor'.", call. = FALSE)
  }
  if (!is.character(sortBy) |
      length(sortBy) != 1 |
      !sortBy %in% c("category", "subcategory", "enrichment", "p_value",
                     "p_adjusted", "q_value", "observed")) {
    stop(paste("'sortBy' must be one of: 'category', 'subcategory',",
               "'enrichment', 'p_value', 'p_adjusted',",
               "'q_value', 'observed'."), call. = FALSE)
  }
  if (!is.logical(convertVer) |
      length(convertVer) != 1) {
    stop("'convertVer' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(mirBaseVer) |
      length(mirBaseVer) != 1) {
    stop(paste("'mirBaseVer' must be a numeric object of lenght one",
               "specifying the miRBase version of miRNA IDs in 'mirnaObj'."),
         call. = FALSE)
  }

  validCategories <- suppressMessages(
    rbioapi::rba_mieaa_cats(mirna_type = mirnaType,
                            species = organism)
  )
  if (!category %in% validCategories & category != "all") {
    stop(paste("'category' must be one of:\n", "\n- all",
               paste(paste("-", validCategories), collapse = "\n"), sep = ""),
         call. = FALSE)
  }

  if (category == "all") {
    cat <- NULL
  } else {
    cat <- category
  }

  ## extract miRNA expression values
  mirna <- mirnaDE(mirnaObj, onlySignificant = FALSE)

  ## order miRNAs on the basis of a specified criterion
  mirnaRes = mirna[order(mirna[, "logFC"], decreasing = TRUE), ]
  mirnaID = mirnaRes$ID

  ## convert miRNAs to miRBase v22
  if (convertVer == TRUE) {
    mirna <- quiet(suppressMessages(
      rbioapi::rba_mieaa_convert_version(mirnaID,
                                         mirna_type = mirnaType,
                                         input_version = mirBaseVer,
                                         output_version = 22,
                                         simple_output = TRUE,
                                         ...)
    ))
  }

  ## perform GSEA on miRNAs
  mirGsea <- quiet(suppressMessages(
    rbioapi::rba_mieaa_enrich(test_set = mirnaID,
                              mirna_type = mirnaType,
                              test_type = "GSEA",
                              species = organism,
                              categories = cat,
                              p_adj_method = pAdjustment,
                              independent_p_adj = independentAdj,
                              sig_level = pCutoff,
                              min_hits = minHits,
                              sort_by = sortBy,
                              ...)
  ))

  ## convert wrong character columns to numeric
  numCols <- c("P-value", "P-adjusted", "Q-value", "Observed")
  if (nrow(mirGsea) > 0) {
    mirGsea[, numCols] <- sapply(mirGsea[, numCols], as.numeric)
  }

  ## rename columns with consistent names
  colnames(mirGsea) <- gsub("-", ".", colnames(mirGsea))

  ## creation of output object
  resGsea <- new("MirnaGsea",
                 result = mirGsea,
                 pvalueCutoff = pCutoff,
                 pAdjustMethod = pAdjustment,
                 organism = organism,
                 ontology = category,
                 gene = mirnaID,
                 lfc = mirnaRes$logFC,
                 keytype = "miRBase v22")
  
  ## report the results of miRNA GSEA
  message(paste("miRNA GSEA reported", nrow(enrichmentResults(resGsea)),
                "significant terms!"))

  ## return GSEA results
  return(resGsea)

}





#' @describeIn MirnaEnrichment-class Subset a
#' [`MirnaEnrichment`][MirnaEnrichment-class] object that contains multiple
#' miEAA 2.0 categories into one that holds just the specified category
#' @export
setMethod("[", c("MirnaEnrichment", "character", "missing", "missing"),
          function(x, i, j, ..., drop = TRUE) {
            res <- enrichmentResults(x)
            enrichmentResults(x) <- res[res$Category == i, ]
            enrichmentDatabase(x) <- i
            x
          })





#' @describeIn MirnaGsea-class Subset a [`MirnaGsea`][MirnaGsea-class] object
#' that contains multiple miEAA 2.0 categories into one that holds just the
#' specified category
#' @export
setMethod("[", c("MirnaGsea", "character", "missing", "missing"),
          function(x, i, j, ..., drop = TRUE) {
            res <- enrichmentResults(x)
            enrichmentResults(x) <- res[res$Category == i, ]
            enrichmentDatabase(x) <- i
            x
          })





#' Get the list of supported organisms for a given database
#'
#' This function provides the list of supported organisms for different
#' databases, namely Gene Ontology (GO), Kyoto Encyclopedia of Genes and
#' Genomes (KEGG), Reactome, DisGeNet and WikiPathways.
#'
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `Reactome`, `DisGeNet` and `WikiPathways`
#'
#' @returns
#' A `character` listing all the supported organisms for the database
#' specified by the user.
#'
#' @examples
#' # get the supported organisms for GO database
#' supportedOrganisms("GO")
#'
#' # get the supported organisms for Reactome
#' supportedOrganisms("Reactome")
#'
#' @references
#' Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the
#' unification of biology. Nat Genet 25, 25–29 (2000).
#' \url{https://doi.org/10.1038/75556}
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' \url{https://doi.org/10.1093/nar/gki072}
#'
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' \url{https://doi.org/10.1093/nar/gkaa1024}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
supportedOrganisms <- function(database) {

  ## check input
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "Reactome", "DisGeNet", "WikiPathways")) {
    stop(paste("'database' must be one of: 'GO' (default), 'KEGG',",
               "'Reactome', 'DisGeNet', 'WikiPathways'."),
         call. = FALSE)
  }

  ## retrieve the list of supported organisms for a given database
  if (database == "GO") {
    db <- "OrgDb"
  } else {
    db <- database
  }
  supp <- paste(convertOrganism(db, "all"), collapse = ", ")

  ## return the supported organisms
  return(paste("Supported organisms for", database, "database are:", supp))

}





#' Perform an enrichment analysis of integrated microRNA targets
#'
#' This function allows to perform an over-representation analysis (ORA) in
#' order to explore the biological effects of microRNA targets that are
#' statistically associated/correlated with DE-miRNAs. The enrichment analysis
#' can be performed using different databases, namely Gene Ontology
#' (GO) - Biological Process, Kyoto Encyclopedia of Genes and Genomes (KEGG),
#' Reactome, DisGeNet and WikiPathways.
#'
#' @details
#' If the enrichment analysis is performed using `GO` database, the
#' [clusterProfiler::simplify()] method contained in the `clusterProfiler`
#' package may be used. This is particularly useful since the structure of GO
#' database is highly redundant. Therefore, to remove some of the redundancy
#' from the resulting enriched GO terms, the parameter `simplifyGO` can be set
#' to `TRUE`.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `Reactome`, `DisGeNet` and `WikiPathways`.
#' Default is `GO`
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param simplifyGO Logical, whether to apply the
#' [clusterProfiler::simplify()] function in `clusterProfiler` package to
#' reduce the redundancy of GO terms in the results (default is `TRUE`).
#' This argument is only considered when `database` is set to `GO`. For further
#' information see the *details* section
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param qCutoff The q-value cutoff to use for statistical significance
#' The default value is `0.2`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param ... Other parameters that can be given to the enrichment functions
#'
#' @returns
#' A `list` object with two elements, namely 'UP-miRNA targets' and
#' 'DOWN-miRNA targets', containing enrichment results of genes targeted by
#' upregulated miRNAs and by downregulated miRNAs, respectively. Each element
#' of this `list` consists of an object of class
#' [`enrichResult-class`][DOSE::enrichResult-class] containing the outcome
#' of the enrichment analysis of target genes. To access the enrichment results
#' as a standard `data.frame`, you can use the `as.data.frame()` function.
#' The user can also refer to the documentation of `clusterProfiler` for
#' analyzing and plotting results:
#' \url{http://yulab-smu.top/biomedical-knowledge-mining-book/index.html}.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform enrichment analysis of integrated targets with Reactome
#' targets_enrichment <- enrichTargets(obj, database = "Reactome")
#'
#' # extract the enrichment of targets of upregulated miRNAs
#' enr_up <- targets_enrichment[["UP-miRNA targets"]]
#'
#' # extract enrichment results as a data.frame
#' enr_df <- as.data.frame(enr_up)
#' enr_df
#'
#' # create a dotplot of enriched terms
#' enrichplot::dotplot(enr_up)
#'
#' @note
#' The over-representation analysis is implemented in this function through
#' the packages `clusterProfiler`, `ReactomePA` and `DOSE`.
#'
#' @references
#' T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X
#' Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool
#' for interpreting omics data. The Innovation. 2021, 2(3):100141
#'
#' Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
#' reactome pathway analysis and visualization. Molecular BioSystems 2016,
#' 12(2):477-479
#'
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the
#' unification of biology. Nat Genet 25, 25–29 (2000).
#' \url{https://doi.org/10.1038/75556}
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' \url{https://doi.org/10.1093/nar/gki072}
#'
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' \url{https://doi.org/10.1093/nar/gkaa1024}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichTargets <- function(mirnaObj,
                          database = "GO",
                          organism = "Homo sapiens",
                          simplifyGO = TRUE,
                          pCutoff = 0.05,
                          qCutoff = 0.2,
                          pAdjustment = "fdr",
                          ...) {

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
  if (max(dim(mirnaTargetsIntegration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'integrateMirnaTargets()'",
               "function. See '?integrateMirnaTargets' for the details."),
         call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "Reactome", "DisGeNet", "WikiPathways")) {
    stop(paste("'database' must be one of: 'GO' (default), 'KEGG',",
               "'Reactome', 'DisGeNet', 'WikiPathways'."),
         call. = FALSE)
  }
  if (database == "GO" & !organism %in% convertOrganism("OrgDb", "all")) {
    stop(paste("For GO database, 'organism' must be one of:",
               paste(convertOrganism("OrgDb", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "KEGG" &
             !organism %in% convertOrganism("KEGG", "all")) {
    stop(paste("For KEGG database, 'organism' must be one of:",
               paste(convertOrganism("KEGG", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "Reactome" &
      !organism %in% convertOrganism("Reactome", "all")) {
    stop(paste("For Reactome database, 'organism' must be one of:",
               paste(convertOrganism("Reactome", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "DisGeNet" & organism != "Homo sapiens") {
    stop("For DisGeNet database, only 'Homo sapiens' is supported",
         call. = FALSE)
  } else if (database == "WikiPathways" &
      !organism %in% convertOrganism("WikiPathways", "all")) {
    stop(paste("For WikiPathways database, 'organism' must be one of:",
               paste(convertOrganism("WikiPathways", "all"), collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(simplifyGO) |
      length(simplifyGO) != 1) {
    stop("'simplifyGO' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(pCutoff) |
      length(pCutoff) != 1 |
      pCutoff > 1 |
      pCutoff < 0) {
    stop("'pCutoff' must be a number between 0 and 1! (default is 0.05)",
         call. = FALSE)
  }
  if (!is.numeric(qCutoff) |
      length(qCutoff) != 1 |
      qCutoff > 1 |
      qCutoff < 0) {
    stop("'qCutoff' must be a number between 0 and 1! (default is 0.2)",
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

  ## set organism id for enrichment analysis
  if (database == "GO") {
    orgDb <- convertOrganism("OrgDb", organism)
  } else {
    orgDb <- convertOrganism("OrgDb", organism)
    organism <- convertOrganism(database, organism)
  }

  ## retrieve targets to enrich
  downMirnaTargets <- selectTargets(mirnaObj,
                                    TRUE,
                                    "downregulated")
  upMirnaTargets <- selectTargets(mirnaObj,
                                  TRUE,
                                  "upregulated")

  ## set the universe
  universe <- rownames(mirnaObj)[["genes"]]

  ## perform the enrichment for targets of up- and downregulated miRNAs
  if (!(length(upMirnaTargets) == 1 & any(upMirnaTargets == ""))) {
    message("Enriching targets of upregulated miRNAs...\n")
    upMirna <- enrichInternal(targets = upMirnaTargets,
                              database = database,
                              simplifyGO = simplifyGO,
                              pCutoff = pCutoff,
                              qCutoff = qCutoff,
                              pAdjustment = pAdjustment,
                              organism = organism,
                              orgDb = orgDb,
                              universe = universe,
                              ...)
  } else {
    upMirna <- NULL
  }

  if (!(length(downMirnaTargets) == 1 & any(downMirnaTargets == ""))) {
    message("Enriching targets of downregulated miRNAs...\n")
    downMirna <- enrichInternal(targets = downMirnaTargets,
                                database = database,
                                simplifyGO = simplifyGO,
                                pCutoff = pCutoff,
                                qCutoff = qCutoff,
                                pAdjustment = pAdjustment,
                                organism = organism,
                                orgDb = orgDb,
                                universe = universe,
                                ...)
  } else {
    downMirna <- NULL
  }

  ## print a message with the results of the enrichment
  rowDown <- ifelse(
    is.null(downMirna),
    0,
    nrow(downMirna@result[downMirna@result$p.adjust < pCutoff &
                            downMirna@result$qvalue < qCutoff, ])
  )

  rowUp <- ifelse(
    is.null(upMirna),
    0,
    nrow(upMirna@result[upMirna@result$p.adjust < pCutoff &
                          upMirna@result$qvalue < qCutoff, ])
  )

  message(paste("\nThe enrichment of targets reported", rowDown,
                "significantly enriched terms for targets of downregulated",
                "miRNAs and", rowUp, "for targets of upregulated miRNAs.\n"))

  ## return a list object with enrichment of up- and downregulated miRNA targets
  enrRes <- list(upMirna, downMirna)
  names(enrRes) <- c("UP-miRNA targets", "DOWN-miRNA targets")
  return(enrRes)

}





#' Perform an over-representation analysis of differentially expressed genes
#'
#' This function allows to enrich upregulated and downregulated genes through
#' an over-representation analysis (ORA). This can be done to investigate the
#' biological effects of gene expression dysregulations. The enrichment
#' analysis can be carried out using different databases, namely Gene Ontology
#' (GO) - Biological Process, Kyoto Encyclopedia of Genes and Genomes (KEGG),
#' Reactome, DisGeNet and WikiPathways.
#'
#' @details
#' If the enrichment analysis is performed using `GO` database,
#' the [clusterProfiler::simplify()] method contained in the `clusterProfiler`
#' package may be used. This is particularly useful since the structure of GO
#' database is highly redundant. Therefore, to remove some of the redundancy
#' from the resulting enriched GO terms, the parameter `simplifyGO` can be set
#' to `TRUE`.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `Reactome`, `DisGeNet` and `WikiPathways`.
#' Default is `GO`
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param simplifyGO Logical, whether to apply the
#' [clusterProfiler::simplify()] function in `clusterProfiler` package to
#' reduce the redundancy of GO terms in the results (default is `TRUE`).
#' This argument is only considered when `database` is set to `GO`. For further
#' information see the *details* section
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param qCutoff The q-value cutoff to use for statistical significance
#' The default value is `0.2`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param ... Other parameters that can be given to the enrichment functions
#'
#' @returns
#' A `list` object with two elements, namely 'upregulated' and 'downregulated',
#' containing enrichment results of upregulated and downregulated genes,
#' respectively. Each element of this `list` consists of an object of class
#' [`enrichResult-class`][DOSE::enrichResult-class] containing the outcome
#' of the enrichment analysis of DE-genes. To access the enrichment results
#' as a standard `data.frame`, you can use the `as.data.frame()` function.
#' The user can also refer to the documentation of `clusterProfiler` for
#' analyzing and plotting results:
#' \url{http://yulab-smu.top/biomedical-knowledge-mining-book/index.html}.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform enrichment analysis of DE genes with Reactome
#' de_enr <- enrichGenes(obj, database = "Reactome")
#'
#' # extract the enrichment of upregulated genes
#' enr_up <- de_enr[["upregulated"]]
#'
#' # extract enrichment results as a data.frame
#' enr_df <- as.data.frame(enr_up)
#' enr_df
#'
#' # create a dotplot of enriched terms
#' enrichplot::dotplot(enr_up)
#'
#' @note
#' The over-representation analysis is implemented in this function through
#' the packages `clusterProfiler`, `ReactomePA` and `DOSE`.
#'
#' @references
#' T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X
#' Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool
#' for interpreting omics data. The Innovation. 2021, 2(3):100141
#'
#' Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
#' reactome pathway analysis and visualization. Molecular BioSystems 2016,
#' 12(2):477-479
#'
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the
#' unification of biology. Nat Genet 25, 25–29 (2000).
#' \url{https://doi.org/10.1038/75556}
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' \url{https://doi.org/10.1093/nar/gki072}
#'
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' \url{https://doi.org/10.1093/nar/gkaa1024}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichGenes <- function(mirnaObj,
                        database = "GO",
                        organism = "Homo sapiens",
                        simplifyGO = TRUE,
                        pCutoff = 0.05,
                        qCutoff = 0.2,
                        pAdjustment = "fdr",
                        ...) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "Reactome", "DisGeNet", "WikiPathways")) {
    stop(paste("'database' must be one of: 'GO' (default), 'KEGG',",
               "'Reactome', 'DisGeNet', 'WikiPathways'."),
         call. = FALSE)
  }
  if (database == "GO" & !organism %in% convertOrganism("OrgDb", "all")) {
    stop(paste("For GO database, 'organism' must be one of:",
               paste(convertOrganism("OrgDb", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "KEGG" &
             !organism %in% convertOrganism("KEGG", "all")) {
    stop(paste("For KEGG database, 'organism' must be one of:",
               paste(convertOrganism("KEGG", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "Reactome" &
             !organism %in% convertOrganism("Reactome", "all")) {
    stop(paste("For Reactome database, 'organism' must be one of:",
               paste(convertOrganism("Reactome", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "DisGeNet" & organism != "Homo sapiens") {
    stop("For DisGeNet database, only 'Homo sapiens' is supported",
         call. = FALSE)
  } else if (database == "WikiPathways" &
             !organism %in% convertOrganism("WikiPathways", "all")) {
    stop(paste("For WikiPathways database, 'organism' must be one of:",
               paste(convertOrganism("WikiPathways", "all"), collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(simplifyGO) |
      length(simplifyGO) != 1) {
    stop("'simplifyGO' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.numeric(pCutoff) |
      length(pCutoff) != 1 |
      pCutoff > 1 |
      pCutoff < 0) {
    stop("'pCutoff' must be a number between 0 and 1! (default is 0.05)",
         call. = FALSE)
  }
  if (!is.numeric(qCutoff) |
      length(qCutoff) != 1 |
      qCutoff > 1 |
      qCutoff < 0) {
    stop("'qCutoff' must be a number between 0 and 1! (default is 0.2)",
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

  ## set organism id for enrichment analysis
  if (database == "GO") {
    orgDb <- convertOrganism("OrgDb", organism)
  } else {
    orgDb <- convertOrganism("OrgDb", organism)
    organism <- convertOrganism(database, organism)
  }

  ## retrieve upregulated and downregulated genes
  de <- geneDE(mirnaObj)
  up <- de$ID[de$logFC > 0]
  down <- de$ID[de$logFC < 0]

  ## set the universe
  universe <- rownames(mirnaObj)[["genes"]]

  ## perform the enrichment for upregulated genes
  if (length(up) >= 1) {
    message("Performing the enrichment of upregulated genes...\n")
    upEnr <- enrichInternal(targets = up,
                            database = database,
                            simplifyGO = simplifyGO,
                            pCutoff = pCutoff,
                            qCutoff = qCutoff,
                            pAdjustment = pAdjustment,
                            organism = organism,
                            orgDb = orgDb,
                            universe = universe,
                            ...)
  } else {
    upEnr <- NULL
  }

  ## perform the enrichment for downregulated genes
  if (length(down) >= 1) {
    message("Performing the enrichment of downregulated genes...\n")
    downEnr <- enrichInternal(targets = down,
                              database = database,
                              simplifyGO = simplifyGO,
                              pCutoff = pCutoff,
                              qCutoff = qCutoff,
                              pAdjustment = pAdjustment,
                              organism = organism,
                              orgDb = orgDb,
                              universe = universe,
                              ...)
  } else {
    downEnr <- NULL
  }

  ## print a message with the results of the enrichment
  rowDown <- ifelse(
    is.null(downEnr),
    0,
    nrow(downEnr@result[downEnr@result$p.adjust < pCutoff &
                          downEnr@result$qvalue < qCutoff, ])
  )

  rowUp <- ifelse(
    is.null(upEnr),
    0,
    nrow(upEnr@result[upEnr@result$p.adjust < pCutoff &
                        upEnr@result$qvalue < qCutoff, ])
  )

  message(paste("\nThe enrichment of targets reported", rowDown,
                "significantly enriched terms for downregulated genes",
                "and", rowUp, "for upregulated genes.\n"))

  ## return a list object with enrichment of up- and downregulated genes
  enrRes <- list(upEnr, downEnr)
  names(enrRes) <- c("upregulated", "downregulated")
  return(enrRes)

}





## enrich genes using clusterProfiler, DOSE and ReactomePA
enrichInternal <- function(targets,
                           database,
                           simplifyGO,
                           pCutoff,
                           qCutoff,
                           pAdjustment,
                           organism,
                           orgDb,
                           universe,
                           ...) {

  ## perform ORA analysis with the selected database
  if (database == "GO") {

    ## enrich targets with GO
    enr <- clusterProfiler::enrichGO(gene = targets,
                                     ont = "BP",
                                     OrgDb = orgDb,
                                     keyType = "SYMBOL",
                                     pAdjustMethod = pAdjustment,
                                     pvalueCutoff = pCutoff,
                                     qvalueCutoff = qCutoff,
                                     universe = universe,
                                     ...)

    ## simplify GO results if wanted
    if (simplifyGO == TRUE) {
      if (!is.null(enr)) enr <- clusterProfiler::simplify(enr)
    }

  } else {

    ## conversion of gene symbols to entrez ID
    targetsEntrez <- clusterProfiler::bitr(geneID = targets,
                                           fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = orgDb)
    univEntrez <- clusterProfiler::bitr(geneID = universe,
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = orgDb)

    ## perform the enrichment with other databases
    if (database == "KEGG") {

      ## enrich targets with KEGG database
      enr <- clusterProfiler::enrichKEGG(gene = targetsEntrez$ENTREZID,
                                         organism = organism,
                                         pAdjustMethod = pAdjustment,
                                         pvalueCutoff = pCutoff,
                                         qvalueCutoff = qCutoff,
                                         universe = univEntrez$ENTREZID,
                                         ...)

    } else if (database == "Reactome") {

      ## enrich targets with Reactome
      enr <- ReactomePA::enrichPathway(gene = targetsEntrez$ENTREZID,
                                       organism = organism,
                                       pAdjustMethod = pAdjustment,
                                       pvalueCutoff = pCutoff,
                                       qvalueCutoff = qCutoff,
                                       universe = univEntrez$ENTREZID,
                                       readable = TRUE,
                                       ...)

    } else if (database == "DisGeNet") {

      ## enrich targets with DisGeNet
      enr <- DOSE::enrichDGN(gene = targetsEntrez$ENTREZID,
                             pAdjustMethod = pAdjustment,
                             pvalueCutoff = pCutoff,
                             qvalueCutoff = qCutoff,
                             universe = univEntrez$ENTREZID,
                             readable = TRUE,
                             ...)

    } else if (database == "WikiPathways") {

      ## enrich targets with Wiki Pathways
      enr <- clusterProfiler::enrichWP(gene = targetsEntrez$ENTREZID,
                                       organism = organism,
                                       pAdjustMethod = pAdjustment,
                                       pvalueCutoff = pCutoff,
                                       qvalueCutoff = qCutoff,
                                       universe = univEntrez$ENTREZID,
                                       ...)
    }
  }

  ## return the enrichment object
  return(enr)

}





#' Perform a gene set enrichment analysis (GSEA)
#'
#' This function allows to carry out a gene set enrichment analysis (GSEA) on
#' gene expression data. This can be done to investigate the biological effects
#' of gene expression dysregulations. The GSEA can be performed using different
#' databases, namely Gene Ontology (GO) - Biological Process, Kyoto
#' Encyclopedia of Genes and Genomes (KEGG), Reactome, DisGeNet
#' and WikiPathways. The ranking parameter used by this function to apply the
#' GSEA algorithm is `logFC`.
#'
#' @details
#' If GSEA is performed using `GO` database, the [clusterProfiler::simplify()]
#' method contained in the `clusterProfiler` package may be used. This is
#' particularly useful since the structure of GO database is highly redundant.
#' Therefore, to remove some of the redundancy from the resulting enriched GO
#' terms, the parameter `simplifyGO` can be set to `TRUE`.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `Reactome`, `DisGeNet` and `WikiPathways`.
#' Default is `GO`
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param simplifyGO Logical, whether to apply the
#' [clusterProfiler::simplify()] function in `clusterProfiler` package to
#' reduce the redundancy of GO terms in the results (default is `TRUE`).
#' This argument is only considered when `database` is set to `GO`. For further
#' information see the *details* section
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param minGSSize The minimum size of genes that a gene set should contain to
#' consider it in the analysis. Default is `10` because smaller gene sets may
#' report statistical significance just with a couple of genes
#' @param maxGSSize The maximum size of genes that a gene set should contain to
#' consider it in the analysis. Default is `500` because bigger gene sets may
#' contain too general categories with a high probability of obtaining
#' statistical significance
#' @param eps The default lower bound for estimating p-values. Default is
#' `1e-10`. To calculate smaller p-values in a more accurate way, set this
#' parameter to `0`
#' @param ... Other parameters that can be given to the GSEA functions
#'
#' @returns
#' An object of class [`gseaResult-class`][DOSE::gseaResult-class]
#' containing the outcome of the gene set enrichment analysis (GSEA). To
#' access the enrichment results as a standard `data.frame`, you can use the
#' `as.data.frame()` function. The user can also refer to the documentation of
#' `clusterProfiler` for analyzing and plotting results:
#' \url{http://yulab-smu.top/biomedical-knowledge-mining-book/index.html}.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform a gene set enrichment analysis with KEGG
#' gse <- gseaGenes(obj, database = "KEGG")
#'
#' # extract GSEA results as a data.frame
#' gse_df <- as.data.frame(gse)
#' gse_df
#'
#' # create a dotplot of enriched terms
#' enrichplot::dotplot(gse)
#'
#' @note
#' The GSEA algorithm is applied in this function through the packages
#' `clusterProfiler`, `ReactomePA` and `DOSE`.
#'
#' @references
#' T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X
#' Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool
#' for interpreting omics data. The Innovation. 2021, 2(3):100141
#'
#' Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
#' reactome pathway analysis and visualization. Molecular BioSystems 2016,
#' 12(2):477-479
#'
#' Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
#' R/Bioconductor package for Disease Ontology Semantic and Enrichment
#' analysis. Bioinformatics 2015 31(4):608-609
#'
#' Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the
#' unification of biology. Nat Genet 25, 25–29 (2000).
#' \url{https://doi.org/10.1038/75556}
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' \url{https://doi.org/10.1093/nar/gki072}
#'
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' \url{https://doi.org/10.1093/nar/gkaa1024}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
gseaGenes <- function(mirnaObj,
                      database = "GO",
                      organism = "Homo sapiens",
                      simplifyGO = TRUE,
                      pCutoff = 0.05,
                      pAdjustment = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 1e-10,
                      ...) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "Reactome", "DisGeNet", "WikiPathways")) {
    stop(paste("'database' must be one of: 'GO' (default), 'KEGG',",
               "'Reactome', 'DisGeNet', 'WikiPathways'."),
         call. = FALSE)
  }
  if (database == "GO" & !organism %in% convertOrganism("OrgDb", "all")) {
    stop(paste("For GO database, 'organism' must be one of:",
               paste(convertOrganism("OrgDb", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "KEGG" &
             !organism %in% convertOrganism("KEGG", "all")) {
    stop(paste("For KEGG database, 'organism' must be one of:",
               paste(convertOrganism("KEGG", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "Reactome" &
             !organism %in% convertOrganism("Reactome", "all")) {
    stop(paste("For Reactome database, 'organism' must be one of:",
               paste(convertOrganism("Reactome", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "DisGeNet" & organism != "Homo sapiens") {
    stop("For DisGeNet database, only 'Homo sapiens' is supported",
         call. = FALSE)
  } else if (database == "WikiPathways" &
             !organism %in% convertOrganism("WikiPathways", "all")) {
    stop(paste("For WikiPathways database, 'organism' must be one of:",
               paste(convertOrganism("WikiPathways", "all"), collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(simplifyGO) |
      length(simplifyGO) != 1) {
    stop("'simplifyGO' must be logical (TRUE/FALSE)!", call. = FALSE)
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
  if (!is.numeric(minGSSize) |
      length(minGSSize) != 1 |
      minGSSize < 0) {
    stop(paste("'minGSSize' must be a number higher than 0",
               "(default is 10). For details see ?gseaGenes"),
         call. = FALSE)
  }
  if (!is.numeric(maxGSSize) |
      length(maxGSSize) != 1 |
      maxGSSize < 0 |
      maxGSSize < minGSSize) {
    stop(paste("'maxGSSize' must be a number higher than 0",
               "and it must be greater than 'minGSSize' (default is 500).",
               "For details see ?gseaGenes"),
         call. = FALSE)
  }
  if (!is.numeric(eps) |
      length(eps) != 1 |
      eps > 1 |
      eps < 0) {
    stop(paste("'eps' must be a very low number that represents the lowest",
               "p-value that can be calculated (default is 1e-10).",
               "For details see ?gseaGenes"),
         call. = FALSE)
  }

  ## set organism id for enrichment analysis
  if (database == "GO") {
    orgDb <- convertOrganism("OrgDb", organism)
  } else {
    orgDb <- convertOrganism("OrgDb", organism)
    organism <- convertOrganism(database, organism)
  }

  ## extract gene differential expression
  de <- geneDE(mirnaObj, onlySignificant = FALSE)

  ## order differential expression results based on logFC
  de <- de[order(de$logFC, decreasing = TRUE), ]

  ## create ranked gene list
  rankGenes <- de$logFC
  names(rankGenes) <- de$ID
  
  ## check if there are NAs in differential expression results
  if (any(is.na(rankGenes))) {
    warning(paste(sum(is.na(rankGenes)), "genes do not have `logFC` values.",
                  "These genes will be excluded from GSEA..."), call. = TRUE)
    rankGenes <- rankGenes[!is.na(rankGenes)]
  }

  ## perform gene set enrichment analysis
  gse <- gseaInternal(genes = rankGenes,
                      database = database,
                      simplifyGO = simplifyGO,
                      pCutoff = pCutoff,
                      pAdjustment = pAdjustment,
                      organism = organism,
                      orgDb = orgDb,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      eps = eps,
                      ...)

  ## return the results of gene set enrichment analysis
  return(gse)

}





# perform GSEA using clusterProfiler, DOSE and ReactomePA
gseaInternal <- function(genes,
                         database,
                         simplifyGO,
                         pCutoff,
                         pAdjustment,
                         organism,
                         orgDb,
                         minGSSize,
                         maxGSSize,
                         eps,
                         ...) {

  ## perform GSEA analysis with the selected database
  if (database == "GO") {

    ## enrich targets with GO
    enr <- clusterProfiler::gseGO(geneList = genes,
                                  ont = "BP",
                                  OrgDb = orgDb,
                                  keyType = "SYMBOL",
                                  minGSSize = minGSSize,
                                  maxGSSize = maxGSSize,
                                  pAdjustMethod = pAdjustment,
                                  pvalueCutoff = pCutoff,
                                  eps = eps,
                                  ...)

    ## simplify GO results if wanted
    if (simplifyGO == TRUE) {
      if (!is.null(enr)) enr <- clusterProfiler::simplify(enr)
    }

  } else {

    ## conversion of gene symbols to entrez ID
    geneEntrez <- clusterProfiler::bitr(geneID = names(genes),
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = orgDb)
    
    ## remove duplicated mappings if present
    geneEntrez <- geneEntrez[which(!duplicated(geneEntrez$SYMBOL)), ]
    
    ## recreate ranked list with entrez ID
    convGenes <- genes[names(genes) %in% geneEntrez$SYMBOL]
    names(convGenes) <- geneEntrez$ENTREZID

    ## perform the enrichment with other databases
    if (database == "KEGG") {

      ## enrich targets with KEGG database
      enr <- clusterProfiler::gseKEGG(geneList = convGenes,
                                      organism = organism,
                                      pAdjustMethod = pAdjustment,
                                      pvalueCutoff = pCutoff,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps = eps,
                                      ...)

    } else if (database == "Reactome") {

      ## enrich targets with Reactome
      enr <- ReactomePA::gsePathway(geneList = convGenes,
                                    organism = organism,
                                    pAdjustMethod = pAdjustment,
                                    pvalueCutoff = pCutoff,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    eps = eps,
                                    ...)

    } else if (database == "DisGeNet") {

      ## enrich targets with DisGeNet
      enr <- DOSE::gseDGN(geneList = convGenes,
                          pAdjustMethod = pAdjustment,
                          pvalueCutoff = pCutoff,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          eps = eps,
                          ...)

    } else if (database == "WikiPathways") {

      ## enrich targets with Wiki Pathways
      enr <- clusterProfiler::gseWP(geneList = convGenes,
                                    organism = organism,
                                    pAdjustMethod = pAdjustment,
                                    pvalueCutoff = pCutoff,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    eps = eps,
                                    ...)

    }
  }

  ## return the enrichment object
  return(enr)

}

