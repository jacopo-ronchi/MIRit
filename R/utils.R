## identify column names
identifyColNames <- function(tabOutput, tabID) {

  # define accepted column names
  idAccepted <- "ID|Symbol|Gene_Symbol|Mirna|mir|Gene|gene.symbol|Gene.symbol"
  fcAccepted <- "logFC|log2FoldChange|FC|lFC"
  exprAccepted <- "AveExpr|baseMean|logCPM"
  pvalAccepted <- "P.Value|pvalue|PValue|Pvalue"
  fdrAccepted <- "adj.P.Val|padj|FDR|fdr|adj|adj.p|adjp"
  acceptedNames <- c(idAccepted, fcAccepted, exprAccepted,
                     pvalAccepted, fdrAccepted)

  ## try to identify column names
  dfNames <- colnames(tabOutput)
  idCol <- grep(idAccepted, dfNames)
  fcCol <- grep(fcAccepted, dfNames)
  exprCol <- grep(exprAccepted, dfNames)
  pvalCol <- grep(pvalAccepted, dfNames)
  fdrCol <- grep(fdrAccepted, dfNames)
  tableCols <- list(idCol, fcCol, exprCol, pvalCol, fdrCol)
  tableNames <- c("ID", "logFC", "AveExpr", "P.Value", "FDR")

  ## check if columns are correctly identified
  if (any(lengths(tableCols) == 0)) {
    stop(paste("The software is unable to automatically find columns",
               "relative to",
               paste(tableNames[which(lengths(tableCols) == 0)],
                     collapse = ", "),
               "in the", tabID, "differential expression data.frame!",
               "Please change the name of these columns to one of:",
               paste(acceptedNames[which(lengths(tableCols) == 0)],
                     collapse = ", ")), call. = FALSE)
  } else if (any(lengths(tableCols) > 1)) {
    stop(paste("More than one column can be interpreted as",
               paste(tableNames[which(lengths(tableCols) > 1)],
                     collapse = ", "),
               "in the", tabID, "differential expression data.frame!",
               "Please rename these columns to have unambiguous column names!",
               "See ?MirnaExperiment for further details"), call. = FALSE)
  }

  ## create a data.frame with desired columns
  tabOutput <- tabOutput[, as.numeric(tableCols)]
  colnames(tabOutput) <- tableNames

  ## return data.frame
  return(tabOutput)

}





## obtain integrated (or not) targets to enrich
selectTargets <- function(mirnaObj, integratedTargets, direction) {

  if (integratedTargets == TRUE) {

    ## check that integration has been performed
    if (nrow(mirnaTargetsIntegration(mirnaObj)) == 0) {
      stop(paste("No targets integrated with DE-miRNAs have been found!",
                 "Run 'integrateMirnaTargets()' to obtain targets",
                 "whose expression is linked to that of DE-miRNAs.",
                 "Instead, if you want to enrich all the targets of DE-miRNAs",
                 "just set 'integratedTargets = FALSE'."), call. = FALSE)
    }

    ## extract integrated targets
    intRes <- mirnaTargetsIntegration(mirnaObj)
    if (colnames(intRes)[4] == "Correlation.Type") {
      targets <- unique(intRes$Target[intRes$microRNA.Direction == direction])
    } else if (colnames(intRes)[4] == "n_NON_DE_targets") {
      targets <- intRes$DE_targets[intRes$direction == direction]
      targets <- paste(targets, collapse = "/")
      targets <- stringr::str_split(targets, "/")
      targets <- unlist(targets)
    }

  } else {

    ## extract differential expression results
    dem <- mirnaDE(mirnaObj)
    if (direction == "upregulated") {
      demId <- dem$ID[dem$logFC > 0]
    } else if (direction == "downregulated") {
      demId <- dem$ID[dem$logFC < 0]
    }

    ## extract all targets of DE-miRNAs
    targets <- mirnaTargets(mirnaObj)
    targets <- targets$target_symbol[targets$mature_mirna_id %in% demId]
    targets <- unique(targets)

  }

  ## return targets
  return(targets)

}





## hide output by cat
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}





## convert organism ids
convertOrganism <- function(col, orgId) {

  ## create a dataframe listing supported organisms for different functions
  supp <- data.frame(multiMiR = c("Homo sapiens", "Mus musculus",
                                  "Rattus norvegicus", NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA),
                     miEAA = c("Homo sapiens", "Mus musculus",
                               "Rattus norvegicus", "Arabidopsis thaliana",
                               "Bos taurus", "Caenorhabditis elegans",
                               "Drosophila melanogaster", "Danio rerio",
                               "Gallus gallus", "Sus scrofa",
                               NA, NA, NA, NA, NA, NA, NA),
                     OrgDb = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
                               "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
                               "org.Dm.eg.db", "org.Dr.eg.db", "org.Gg.eg.db",
                               "org.Ss.eg.db", "org.Sc.sgd.db", NA,
                               "org.Pt.eg.db", NA, NA, "org.Cf.eg.db",
                               "org.Ag.eg.db"),
                     KEGG = c("hsa", "mmu", "rno", "ath", "bta", "cel", "dme",
                              "dre", "gga", "ssc", "sce", "pop", "ptr", "osa",
                              "ecb", "cfa", "aga"),
                     Reactome = c("human", "mouse", "rat", NA, NA, "celegans",
                                  "fly", "zebrafish", NA, NA, "yeast", NA, NA,
                                  NA, NA, NA, NA),
                     DisGeNet = c("Homo sapiens", NA, NA, NA, NA, NA, NA, NA,
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA),
                     WikiPathways = c("Homo sapiens", "Mus musculus",
                                      "Rattus norvegicus",
                                      "Arabidopsis thaliana",
                                      "Bos taurus", "Caenorhabditis elegans",
                                      "Drosophila melanogaster",
                                      "Danio rerio", "Gallus gallus",
                                      "Sus scrofa", "Saccharomyces cerevisiae",
                                      "Populus trichocarpa", "Pan troglodytes",
                                      "Oryza sativa", "Equus caballus",
                                      "Canis familiaris", "Anopheles gambiae"),
                     graph_kegg = c("hsapiens", "mmusculus", "rnorvegicus",
                                    "athaliana", "btaurus", "celegans",
                                    "dmelanogaster", "drerio", "ggallus",
                                    "sscrofa", "scerevisiae", NA, NA, NA, NA,
                                    "cfamiliaris", NA),
                     graph_reactome = c("hsapiens", "mmusculus", "rnorvegicus",
                                        NA, "btaurus", "celegans",
                                        "dmelanogaster", "drerio", "ggallus",
                                        "sscrofa", "scerevisiae", NA, NA, NA,
                                        NA, "cfamiliaris", NA),
                     graph_wikipathways = c("hsapiens", "mmusculus",
                                            "rnorvegicus", "athaliana",
                                            "btaurus", "celegans",
                                            "dmelanogaster", "drerio",
                                            "ggallus", "sscrofa", "scerevisiae",
                                            NA, NA, NA, NA, "cfamiliaris", NA))
  rownames(supp) <- supp$WikiPathways

  ## convert organism ID or return all supported organisms
  if (orgId != "all") {
    id <- supp[orgId, col]
  } else {
    id <- rownames(supp[which(!is.na(supp[, col])), ])
  }

  ## return id
  return(id)

}

