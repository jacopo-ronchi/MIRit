## helper function to download the desired gene set for enrichment analyses
prepareGeneSet <- function(organism, database, category) {
  
  ## determine organism name accepted by database
  org <- species[species$specie == organism, database]
  
  ## determine organism name for ID conversion
  convOrg <- species[species$specie == organism, "Conversion"]
  
  ## download and prepare the appropriate gene set
  if (database == "GO") {
    gs <- geneset::getGO(org = org,
                         ont = category)
    colnames(gs$geneset) <- c("id", "symbol")
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
  } else if (database == "KEGG") {
    gs <- geneset::getKEGG(org = org,
                           category = category)
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
    symb <- genekitr::transId(gs$gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
  } else if (database == "MsigDB") {
    gs <- geneset::getMsigdb(org = org,
                             category = category)
    gs <- gs$geneset
    symb <- genekitr::transId(gs$entrez_gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$entrez_gene, symb$input_id)]
    colnames(gs)[1] <- "name"
  } else if (database == "WikiPathways") {
    gs <- geneset::getWiki(org = org)
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
    symb <- genekitr::transId(gs$gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
  } else if (database == "Reactome") {
    gs <- geneset::getReactome(org = org)
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
    symb <- genekitr::transId(gs$gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
  } else if (database == "Enrichr") {
    gs <- geneset::getEnrichrdb(org = org,
                                library = category)
    gs <- gs$geneset
    colnames(gs) <- c("name", "symbol")
  } else if (database == "DO") {
    gs <- geneset::getHgDisease(source = "do")
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
    symb <- genekitr::transId(gs$gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
    dupDO <- duplicated(data.frame("name" = tolower(gs$name),
                                   "symbol" = gs$symbol))
    gs <- gs[!dupDO, ]
  } else if (database == "NCG") {
    if (category == "v6") {
      gs <- geneset::getHgDisease(source = "ncg_v6")
      gs <- gs$geneset
      symb <- genekitr::transId(gs$gene,
                                transTo = "symbol",
                                org = convOrg)
      gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
      colnames(gs)[1] <- "name"
    } else if (category == "v7") {
      gs <- geneset::getHgDisease(source = "ncg_v7")
      gs <- gs$geneset
      symb <- genekitr::transId(gs$gene,
                                transTo = "symbol",
                                org = convOrg)
      gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
      colnames(gs)[1] <- "name"
    }
  } else if (database == "DisGeNET") {
    gs <- geneset::getHgDisease(source = "disgenet")
    gs <- merge(gs$geneset, gs$geneset_name, by = "id")
    symb <- genekitr::transId(gs$gene,
                              transTo = "symbol",
                              org = convOrg)
    gs$symbol <- symb$symbol[match(gs$gene, symb$input_id)]
  } else if (database == "COVID19") {
    gs <- geneset::getHgDisease(source = "covid19")
    gs <- gs$geneset
    colnames(gs) <- c("name", "symbol")
  }
  
  ## convert gene set to list
  gs <- split(gs$symbol, gs$name)
  
  ## return gene set
  return(gs)
  
}





## helper function that checks for valid categories for a given database
validateCategories <- function(database, category, organism) {
  
  ## check if category is included in the specified database
  if (database == "GO" &
      !category %in% c("bp", "mf", "cc")) {
    stop(paste("For GO database, 'category' must be one of 'bp', 'mf', 'cc'.",
               "For additional details, see ?enrichGenes"))
  } else if (database == "KEGG" &
             !category %in% c("pathway", "module", "enzyme",
                              "disease", "drug", "network")) {
    stop(paste("For KEGG database, 'category' must be one of 'pathway',",
               "'module', 'enzyme', 'disease', 'drug', 'network'.",
               "For additional details, see ?enrichGenes"))
  } else if(database == "KEGG" &
            category %in% c("disease", "drug", "network") &
            organism != "Homo sapiens") {
    stop(paste("For KEGG database, the categories 'disease', 'drug', and",
               "'network' are available only for specie Homo sapiens.",
               "For additional details, see ?enrichGenes"))
  } else if (database == "MsigDB" &
             !category %in% c("H", "C1", "C2-CGP", "C2-CP-BIOCARTA",
                              "C2-CP-KEGG", "C2-CP-PID", "C2-CP-REACTOME",
                              "C2-CP-WIKIPATHWAYS", "C3-MIR-MIRDB",
                              "C3-MIR-MIR_Legacy", "C3-TFT-GTRD",
                              "C3-TFT-TFT_Legacy", "C4-CGN", "C4-CM",
                              "C5-GO-BP", "C5-GO-CC", "C5-GO-MF", "C5-HPO",
                              "C6", "C7-IMMUNESIGDB", "C7-VAX", "C8")) {
    stop(paste("For MsigDB database, 'category' must be one of 'H', 'C1',",
               "'C2-CGP', 'C2-CP-BIOCARTA', 'C2-CP-KEGG', 'C2-CP-PID',",
               "'C2-CP-REACTOME', 'C2-CP-WIKIPATHWAYS', 'C3-MIR-MIRDB'",
               "'C3-MIR-MIR_Legacy', 'C3-TFT-GTRD', 'C3-TFT-TFT_Legacy',",
               "'C4-CGN', 'C4-CM', 'C5-GO-BP', 'C5-GO-CC', 'C5-GO-MF',",
               "'C5-HPO', 'C6', 'C7-IMMUNESIGDB', 'C7-VAX', 'C8'.",
               "For additional details, see ?enrichGenes"))
  } else if (database == "NCG" &
             !category %in% c("v6", "v7")) {
    stop(paste("For NCG database, 'category' must be one of 'v6', 'v7'.",
               "For additional details, see ?enrichGenes"))
  } else if (database == "Enrichr" &
             !category %in% geneset::enrichr_metadata$library[
               geneset::enrichr_metadata$organism == species$Enrichr[
                 species$specie == organism]]) {
    stop(paste("Valid categories for Enrichr database are listed in",
               "'geneset::enrichr_metadata'. For additional details,",
               "see ?enrichGenes"))
    
  }
}





#' Perform functional enrichment analysis of genes
#'
#' This function allows to investigate the biological functions and pathways
#' that result dysregulated across biological conditions. In particular,
#' different enrichment approaches can be used, including over-representation
#' analysis (ORA), gene-set enrichment analysis (GSEA), and Correlation Adjusted
#' MEan RAnk gene set test (CAMERA). Moreover, for all these analyses, the
#' enrichment can be carried out using different databases, namely Gene Ontology
#' (GO), Kyoto Encyclopedia of Genes and Genomes (KEGG), MsigDB, WikiPathways,
#' Reactome, Enrichr, Disease Ontology (DO), Network of Cancer Genes (NCG),
#' DisGeNET, and COVID19. For exhaustive information on how to use this
#' function, please refer to the *details* section.
#'
#' @details
#' 
#' ## Enrichment method
#' 
#' The method used for functional enrichment analysis will drastically
#' influence the biological results, and thus, it must be carefully chosen.
#' `ORA` (Boyle et al., 2004) takes differentially expressed genes (separately
#' considering upregulated and downregulated features) and uses the
#' hypergeometric test to infer the biological processes that are regulated by
#' these genes more than would be expected by chance. The downside of this
#' approach is that we only consider genes that passed a pre-defined threshold,
#' thus losing all the slight changes in gene expression that may have important
#' biological consequences.
#' 
#' To address this limit, `GSEA` was introduced (Subramanian, 2005). This
#' analysis starts by ranking genes according to a specific criterion, and then
#' uses a running statistic that is able to identify even slight but coordinated
#' expression changes of genes belonging to a specific pathway. Therefore,
#' `GSEA` is the default method used in MIRit to perform the functional
#' enrichment analysis of genes.
#' 
#' Moreover, in addition to ORA and GSEA, this function allows to perform the
#' enrichment analysis through `CAMERA` (Wu and Smyth, 2012), which is another
#' competitive test used for functional enrichment of genes. The main advantage
#' of this method is that it adjusts the gene set test statistic according to
#' inter-gene correlations. This is particularly interesting since it was
#' demonstrated that inter-gene correlations may affect the reliability of
#' functional enrichment analyses.
#' 
#' ## Databases and categories
#' 
#' Regarding gene sets, multiple databases can be used to investigate the
#' consequences of gene expression alterations. However, different databases
#' also includes several subcategories with different annotations. To
#' specifically query desired categories, the `category` parameter is used. As
#' a reference, here are listed the available categories for the different
#' databases supported:
#' 
#' * Gene Ontology (GO):
#' 
#'    + `bp`, for GO - Biological Processes;
#'    + `mf`, for GO - Molecular Function;
#'    + `cc`, for GO - Cellular Component;
#' 
#' * Kyoto Encyclopedia of Genes and Genomes (KEGG):
#' 
#'    + `pathway`, for KEGG biological pathways;
#'    + `module`, for KEGG reaction modules;
#'    + `enzyme`, for KEGG enzyme nomenclature;
#'    + `disease`, for KEGG diseases (only Homo sapiens supported);
#'    + `drug`, for KEGG drug targets (only Homo sapiens supported);
#'    + `network`, for KEGG disease/drug perturbation netowrks (only Homo
#'    sapiens supported);
#' 
#' * MsigDB:
#' 
#'    + `H`, for MsigDB hallmark genes of specific biological states/processes;
#'    + `C1`, for gene sets of human chromosome cytogenetic bands;
#'    + `C2-CGP`, for expression signatures of genetic and chemical
#'    perturbations;
#'    + `C2-CP-BIOCARTA`, for canonical pathways gene sets derived from the
#'    BioCarta pathway database;
#'    + `C2-CP-KEGG`, for canonical pathways gene sets derived from the
#'    KEGG pathway database;
#'    + `C2-CP-PID`, for canonical pathways gene sets derived from the
#'    PID pathway database;
#'    + `C2-CP-REACTOME`, for canonical pathways gene sets derived from the
#'    Reactome pathway database;
#'    + `C2-CP-WIKIPATHWAYS`, for canonical pathways gene sets derived from the
#'    WikiPathways database;
#'    + `C3-MIR-MIRDB`, for gene sets containing high-confidence gene-level
#'    predictions of human miRNA targets as catalogued by miRDB v6.0 algorithm; 
#'    + `C3-MIR-MIR_Legacy`, for older gene sets that contain genes sharing
#'    putative target sites of human mature miRNA in their 3'-UTRs;
#'    + `C3-TFT-GTRD`, for genes that share GTRD predicted transcription factor
#'    binding sites in the region -1000,+100 bp around the TSS for the
#'    indicated transcription factor;
#'    + `C3-TFT-TFT_Legacy`, for older gene sets that share upstream
#'    cis-regulatory motifs which can function as potential transcription factor
#'    binding sites;
#'    + `C4-CGN`, for gene sets defined by expression neighborhoods centered
#'    on 380 cancer-associated genes;
#'    + `C4-CM`, for cancer modules as defined by Segal et al. 2004;
#'    + `C5-GO-BP`, for GO - biological process ontology;
#'    + `C5-GO-CC`, for GO - cellular component ontology;
#'    + `C5-GO-MF`, for GO - molecular function ontology;
#'    + `C5-HPO`, for Human Phenotype ontology (HPO);
#'    + `C6`, for gene sets that represent signatures of cellular pathways which
#'    are often dis-regulated in cancer;
#'    + `C7-IMMUNESIGDB`, for manually curated gene sets representing chemical
#'    and genetic perturbations of the immune system;
#'    + `C7-VAX`, for gene sets deriving from the Human Immunology Project
#'    Consortium (HIPC) describing human transcriptomic immune responses to
#'    vaccinations;
#'    + `C8`, for gene sets that contain curated cluster markers for cell types;
#'    
#' * WikiPathways;
#' 
#' * Reactome;
#' 
#' * Enrichr:
#' 
#'    + All avaliable gene sets can be listed through
#'    `geneset::enrichr_metadata`
#' 
#' * Disease Ontology (DO);
#' 
#' * Network of Cancer Genes (NCG):
#' 
#'    + `v6`, for the sixth version;
#'    + `v7`, for the seventh version;
#' 
#' * DisGeNET;
#' 
#' * COVID-19.
#' 
#' ## Supported organisms
#' 
#' For each database, different organisms are supported. To check the supported
#' organisms for a given database, MIRit provides the [supportedOrganism()]
#' function.
#' 
#' ## GSEA ranking statistic
#' 
#' The ranking statistic used to order genes before conducting GSEA is able
#' to influence the biological interpretation of functional enrichment
#' results. Several metrics have been used in scientific literature. MIRit
#' implements the possibility of using `signed.pval`, `logFC`, and `log.pval`.
#' In particular, the simplest option is to rank genes according to their
#' `logFC` value. However, this procedure is biased by higher variance for
#' lowly abundant genes.
#' 
#' Therefore, we recommend to use the `signed.pval` metric, which consists in
#' the p-value of a gene multiplied for the sign of its logFC, i.e.
#' `sign(logFC) * p-value`. Alternatively, `log,pval` metric, which consist in
#' the product of logFC and p-value, i.e. `logFC * p-value` can also be used.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param method The functional enrichment analysis to perform. It must be one
#' of `ORA`, `GSEA` (default), and `CAMERA`. For additional information, see
#' the *details* section
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `MsigDB`, `WikiPathways`, `Reactome`,
#' `Enrichr`, `DO`, `NCG`, `DisGeNET`, `COVID19`. Default is `GO`
#' @param category The desired subcategory of gene sets present in `database`.
#' Please, see the *details* section to check the available categories for
#' each database. Default is NULL to use default categories
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param minSize The minimum size for a gene set. All gene sets containing
#' less than this number of genes will not be considered. Default is 10
#' @param maxSize The maximum size for a gene set. All gene sets containing
#' more than this number of genes will not be considered. Default is 500
#' @param rankMetric The ranking statistic used to order genes before performing
#' GSEA. It must be one of `signed.pval` (default), `logFC`, and `log.pval`.
#' For additional information, refer to the *details* section
#' @param eps The lower boundary for p-value calculation (default is 1e-50).
#' To compute exact p-values, this parameter can be set to 0, even though
#' the analysis will be slower
#'
#' @returns
#' For method `GSEA` and `CAMERA`, this function produces an object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing enrichment
#' results. Instead, when `ORA` is used, this function returns a `list` object
#' with two elements, namely 'upregulated' and 'downregulated',
#' each containing a [`FunctionalEnrichment`][FunctionalEnrichment-class]
#' object storing enrichment results of upregulated and downregulated genes,
#' respectively.
#' 
#' To access results of [`FunctionalEnrichment`][FunctionalEnrichment-class]
#' objects, the user can use the [enrichmentResults()] function. Additionally,
#' MIRit provides several functions to graphically represent enrichment
#' analyses, including [enrichmentBarplot()], [enrichmentDotplot()],
#' [gseaPlot()], and [gseaRidgeplot()].
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG
#' de_enr <- enrichGenes(obj, database = "KEGG")
#'
#' # extract results
#' de_df <- enrichmentResults(de_enr)
#'
#' # create a dotplot of enriched terms
#' enrichmentDotplot(de_enr)
#'
#' @note
#' To download gene sets from the above mentioned databases, MIRit uses the
#' `geneset` R package. Moreover, to perform ORA and GSEA, MIRit implements the
#' `fgsea` algorithm, whereas for CAMERA, the `limma` package is used.
#'
#' @references
#' Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr R
#' package and web server. BMC Bioinformatics 24, 214 (2023).
#' \url{https://doi.org/10.1186/s12859-023-05342-9}.
#' 
#' Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment
#' analysis.” bioRxiv. doi:10.1101/060012,
#' \url{http://biorxiv.org/content/early/2016/06/20/060012}.
#' 
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies.” Nucleic Acids Research, 43(7), e47. \url{doi:10.1093/nar/gkv007}.
#' 
#' Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene
#' correlation. Nucleic Acids Res. 2012 Sep 1;40(17):e133.
#' doi: \url{10.1093/nar/gks461}. Epub 2012 May 25. PMID: 22638577; PMCID:
#' PMC3458527.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichGenes <- function(mirnaObj,
                        method = "GSEA",
                        database = "GO",
                        category = NULL,
                        organism = "Homo sapiens",
                        pCutoff = 0.05,
                        pAdjustment = "fdr",
                        minSize = 10,
                        maxSize = 500,
                        rankMetric = "signed.pval",
                        eps = 1e-50) {
  
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
      !database %in% c("GO", "KEGG", "MsigDB", "WikiPathways", "Reactome",
                       "Enrichr", "DO", "NCG", "DisGeNET", "COVID19")) {
    stop(paste("'database' must be one of 'GO', 'KEGG', 'MsigDB',",
               "'WikiPathways', 'Reactome', 'Enrichr', 'DO', 'NCG',",
               "'DisGeNET', 'COVID19'. For additional details,",
               "see ?enrichGenes"),
         call. = FALSE)
  }
  if (!is.character(method) |
      length(method) != 1 |
      !method %in% c("ORA", "GSEA", "CAMERA")) {
    stop(paste("'method' must be one of 'ORA', 'GSEA', 'CAMERA'.",
               "For additional details, see ?enrichGenes"),
         call. = FALSE)
  }
  if (geneDE(mirnaObj, param = TRUE)$method == "Manually added" &
      method == "CAMERA") {
    stop(paste("Functional enrichment analysis with CAMERA is not available",
               "for user-supplied differential expression results..."),
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
  if (!is.numeric(minSize) |
      length(minSize) != 1 |
      minSize < 0 |
      !is.numeric(maxSize) |
      length(maxSize) != 1 |
      maxSize < 0 |
      minSize > maxSize) {
    stop(paste("'minSize' and 'maxSize' must be two positive numbers!",
               "Additionally, maxSize should be bigger than minSize.",
               "For additional details, see ?enrichGenes"),
         call. = FALSE)
  }
  if (!is.character(rankMetric) |
      length(rankMetric) != 1 |
      !rankMetric %in% c("signed.pval", "logFC", "log.pval")) {
    stop(paste("'rankMetric' must be one of 'signed.pval', 'logFC',",
               "'log.pval'. For additional details, see ?enrichGenes"),
         call. = FALSE)
  }
  if (!is.numeric(eps) |
      length(eps) != 1 |
      eps > 1 |
      eps < 0) {
    stop("'eps' must be a number between 0 and 1! (default is 1e-50)",
         call. = FALSE)
  }
  
  ## check if database is supported for the given specie
  supp <- species[!is.na(species[, database]), "specie"]
  if (!organism %in% supp) {
    stop(paste("For", database, "database, 'organism' must be one of:",
               paste(supp, collapse = ", ")))
  }
  
  ## if NULL, set default categories
  if (is.null(category)) {
    if (!database %in% c("Reactome", "WikiPathways", "DO",
                         "DisGeNET", "COVID19")) {
      if (database == "GO") {
        category <- "bp"
      } else if (database == "KEGG") {
        category <- "pathway"
      } else if (database == "MsigDB") {
        category <- "H"
      } else if (database == "NCG") {
        category <- "v7"
      }
      message(paste("Since not specified, 'category' for", database,
                    "database is set to", category, "(default)."))
    } else {
      category <- ""
    }
  }
  
  ## check if category is included in the specified database
  validateCategories(database, category, organism)
  
  ## download the appropriate gene set
  message("Preparing the appropriate gene set...")
  gs <- prepareGeneSet(organism = organism,
                       database = database,
                       category = category)
  
  ## set a describer for the database used
  if (database %in% c("Reactome", "WikiPathways", "DO",
                      "DisGeNET", "COVID19")) {
    dataInfo <- database
  } else {
    dataInfo <- paste(database, " (category: ", category, ")", sep = "")
  }
  
  ## perform the desired functional enrichment analysis
  if (method == "ORA") {
    res <- oraInternal(mirnaObj = mirnaObj,
                       geneSet = gs,
                       pCutoff = pCutoff,
                       pAdjustment = pAdjustment,
                       minSize = minSize,
                       maxSize = maxSize,
                       dataInfo = dataInfo,
                       organism = organism)
  } else if (method == "GSEA") {
    res <- gseaInternal(mirnaObj = mirnaObj,
                        geneSet = gs,
                        pCutoff = pCutoff,
                        pAdjustment = pAdjustment,
                        minSize = minSize,
                        maxSize = maxSize,
                        dataInfo = dataInfo,
                        organism = organism,
                        rankMetric = rankMetric,
                        eps = eps)
  } else if (method == "CAMERA") {
    res <- cameraInternal(mirnaObj = mirnaObj,
                          geneSet = gs,
                          pCutoff = pCutoff,
                          pAdjustment = pAdjustment,
                          minSize = minSize,
                          maxSize = maxSize,
                          dataInfo = dataInfo,
                          organism = organism)
  }
  
  ## return the results of functional enrichment
  return(res)
  
}





## perform a simple over-representation analysis (ORA)
oraInternal <- function(mirnaObj,
                        geneSet,
                        pCutoff,
                        pAdjustment,
                        minSize,
                        maxSize,
                        dataInfo,
                        organism,
                        integrated = FALSE) {
  
  ## retrieve upregulated and downregulated genes/integrated targets
  if (integrated == TRUE) {
    up <- selectTargets(mirnaObj, miRNA.Direction = "downregulated")
    down <- selectTargets(mirnaObj, miRNA.Direction = "upregulated")
  } else {
    de <- geneDE(mirnaObj)
    up <- de$ID[de$logFC > 0]
    down <- de$ID[de$logFC < 0]
  }
  
  ## set the universe
  universe <- rownames(mirnaObj)[["genes"]]
  
  ## perform over-representation analysis for upregulated genes
  message("Performing the enrichment of upregulated genes...")
  oraUp <- fgsea::fora(pathways = geneSet,
                       genes = up,
                       universe = universe,
                       minSize = minSize,
                       maxSize = maxSize)
  
  ## perform over-representation analysis for downregulated genes
  message("Performing the enrichment of downregulated genes...")
  oraDown <- fgsea::fora(pathways = geneSet,
                         genes = down,
                         universe = universe,
                         minSize = minSize,
                         maxSize = maxSize)
  
  ## adjust p-values through the desired approach
  oraUp$padj <- stats::p.adjust(oraUp$pval, method = pAdjustment)
  oraDown$padj <- stats::p.adjust(oraDown$pval, method = pAdjustment)
  
  ## restrict the results to significant terms
  oraUp <- oraUp[oraUp$padj < pCutoff, ]
  oraDown <- oraDown[oraDown$padj < pCutoff, ]
  
  ## order results based on padj
  oraUp <- oraUp[order(oraUp$padj), ]
  oraDown <- oraDown[order(oraDown$padj), ]
  
  ## inform the user about ORA results
  message(paste("The enrichment of genes reported", nrow(oraDown),
                "significantly enriched terms for downregulated genes",
                "and", nrow(oraUp), "for upregulated genes."))
  
  ## store results as FunctionalEnrichment objects
  upEnr <- new("FunctionalEnrichment",
               data = oraUp,
               method = "Over-Representation Analysis (ORA)",
               organism = organism,
               database = dataInfo,
               pCutoff = pCutoff,
               pAdjustment = pAdjustment,
               features = up,
               statistic = numeric(),
               universe = universe,
               geneSet = geneSet)
  downEnr <- new("FunctionalEnrichment",
                 data = oraDown,
                 method = "Over-Representation Analysis (ORA)",
                 organism = organism,
                 database = dataInfo,
                 pCutoff = pCutoff,
                 pAdjustment = pAdjustment,
                 features = down,
                 statistic = numeric(),
                 universe = universe,
                 geneSet = geneSet)
  
  ## return a list object with both objects
  res <- list(upregulated = upEnr,
              downregulated = downEnr)
  return(res)
  
}





## perform a gene-set enrichment analysis (GSEA)
gseaInternal <- function(mirnaObj,
                         geneSet,
                         pCutoff,
                         pAdjustment,
                         minSize,
                         maxSize,
                         dataInfo,
                         organism,
                         rankMetric,
                         eps) {
  
  ## retrieve gene differential expression
  de <- geneDE(mirnaObj, onlySignificant = FALSE)
  
  ## create GSEA ranked list according to the chosen metric
  message(paste("Ranking genes based on ", rankMetric, "...", sep = ""))
  if (rankMetric == "signed.pval") {
    de$stat <- sign(de$logFC) * (-log10(de$P.Value))
    de <- de[order(de$stat, decreasing = TRUE), ]
    stat <- de$stat
    names(stat) <- de$ID
  } else if (rankMetric == "logFC") {
    de <- de[order(de$logFC, decreasing = TRUE), ]
    stat <- de$logFC
    names(stat) <- de$ID
  } else if (rankMetric == "log.pval") {
    de$stat <- de$logFC * (-log10(de$P.Value))
    de <- de[order(de$stat, decreasing = TRUE), ]
    stat <- de$stat
    names(stat) <- de$ID
  }
  
  ## perform GSEA
  message("Performing gene-set enrichment analysis (GSEA)...")
  gse <- fgsea::fgsea(pathways = geneSet,
                      stats = stat,
                      minSize = minSize,
                      maxSize = maxSize,
                      eps = eps)
  
  ## adjust p-values through the desired approach
  gse$padj <- stats::p.adjust(gse$pval, method = pAdjustment)
  
  ## restrict the results to significant terms
  gse <- gse[gse$padj < pCutoff, ]
  
  ## order results based on padj
  gse <- gse[order(gse$padj), ]
  
  ## inform the user about GSEA results
  message(paste("GSEA reported", nrow(gse), "significantly enriched terms."))
  
  ## store results as a FunctionalEnrichment object
  gseObj <- new("FunctionalEnrichment",
                data = gse,
                method = "Gene-Set Enrichment Analysis (GSEA)",
                organism = organism,
                database = dataInfo,
                pCutoff = pCutoff,
                pAdjustment = pAdjustment,
                features = names(stat),
                statistic = stat,
                universe = character(),
                geneSet = geneSet)
  
  ## return the object
  return(gseObj)
  
}





## perform a competitive gene set test accounting for inter-gene
## correlation (CAMERA)
cameraInternal <- function(mirnaObj,
                           geneSet,
                           pCutoff,
                           pAdjustment,
                           minSize,
                           maxSize,
                           dataInfo,
                           organism) {
  
  ## determine gene set sizes
  sizes <- unlist(lapply(geneSet, length))
  
  ## remove categories that are lowly or overly represented
  geneSet <- geneSet[sizes > minSize &
                       sizes < maxSize]
  
  ## extract gene differential expression results
  de <- geneDE(mirnaObj, param = TRUE)
  
  ## access sample metadata
  meta <- MultiAssayExperiment::colData(mirnaObj)
  meta <- meta[!is.na(meta$geneCol), ]
  
  ## determine the appropriate expression matrix and the experimental design
  if (de$method == "limma") {
    expr <- mirnaObj[["genes"]]
    des <- stats::model.matrix(de$design, data = meta)
  } else if (de$method == "edgeR" |
             de$method == "limma-voom") {
    expr <- geneDE(mirnaObj, returnObject = TRUE)
    des <- expr$design
  } else if (de$method == "DESeq2") {
    message("Applying 'limma-voom' pipeline before using CAMERA...")
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
  message("Performing Correlation Adjusted MEan RAnk gene set test (CAMERA)...")
  if (de$method == "edgeR") {
    
    ## perform miRNA-target integration through 'CAMERA'
    rs <- edgeR::camera.DGEList(y = expr,
                                index = geneSet,
                                design = des,
                                contrast = con)
    
  } else {
    
    ## perform miRNA-target integration through 'CAMERA'
    rs <- limma::camera(y = expr,
                        index = geneSet,
                        design = des,
                        contrast = con)
    
  }
  
  ## adjust p-values according to the specified method
  rs$FDR <- stats::p.adjust(rs$PValue, method = pAdjustment)
  
  ## reshape the resulting data.frame
  rs$pathway <- rownames(rs)
  rs$size <- sizes[match(rs$pathway, names(sizes))]
  rs <- rs[, c(5, 2, 1, 6, 3, 4)]
  colnames(rs) <- c("pathway", "direction", "overlap", "size", "pval", "padj")
  
  ## restrict the results to significant terms
  rs <- na.omit(rs)
  rs <- rs[rs$padj < pCutoff, ]
  
  ## order results based on padj
  rs <- rs[order(rs$padj), ]
  
  ## inform the user about GSEA results
  message(paste("CAMERA reported", nrow(rs), "significantly enriched terms."))
  
  ## store results as a FunctionalEnrichment object
  rsObj <- new("FunctionalEnrichment",
               data = rs,
               method = "Correlation Adjusted MEan RAnk gene set test (CAMERA)",
               organism = organism,
               database = dataInfo,
               pCutoff = pCutoff,
               pAdjustment = pAdjustment,
               features = character(),
               statistic = numeric(),
               universe = character(),
               geneSet = geneSet)
  
  ## return the object
  return(rsObj)
  
}





#' Perform an enrichment analysis of integrated microRNA targets
#'
#' This function allows to perform an over-representation analysis (ORA) of
#' integrated miRNA targets in order to explore the biological effects of
#' targets that are statistically associated/correlated with DE-miRNAs. The
#' enrichment analysis can be performed using different databases, namely
#' Gene Ontology (GO), Kyoto Encyclopedia of Genes and Genomes (KEGG), MsigDB,
#' WikiPathways, Reactome, Enrichr, Disease Ontology (DO), Network of Cancer
#' Genes (NCG), DisGeNET, and COVID19.
#'
#' @details
#' For each database, different organisms are supported. To check the supported
#' organisms for a given database, MIRit provides the [supportedOrganism()]
#' function.
#' 
#' Moreover, since different database support multiple subcategories, the
#' `category` parameter can be set to specify the desired resource. For
#' specific information regarding the available categories for the different
#' databases, check the *details* section of the [enrichGenes()] documentation.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database used for the enrichment analysis.
#' It must be one of: `GO`, `KEGG`, `MsigDB`, `WikiPathways`, `Reactome`,
#' `Enrichr`, `DO`, `NCG`, `DisGeNET`, `COVID19`. Default is `GO`
#' @param category The desired subcategory of gene sets present in `database`.
#' Please, see the *details* section to check the available categories for
#' each database. Default is NULL to use default categories
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param minSize The minimum size for a gene set. All gene sets containing
#' less than this number of genes will not be considered. Default is 10
#' @param maxSize The maximum size for a gene set. All gene sets containing
#' more than this number of genes will not be considered. Default is 500
#'
#' @returns
#' This function produces a `list` object with two elements, namely
#' 'upregulated' and 'downregulated', each containing a
#' [`FunctionalEnrichment`][FunctionalEnrichment-class]
#' object storing enrichment results of upregulated and downregulated target
#' genes, respectively.
#' 
#' To access results of [`FunctionalEnrichment`][FunctionalEnrichment-class]
#' objects, the user can use the [enrichmentResults()] function. Additionally,
#' MIRit provides several functions to graphically represent enrichment
#' analyses, including [enrichmentBarplot()], and [enrichmentDotplot()].
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform enrichment analysis of integrated targets with KEGG
#' targets_enrichment <- enrichTargets(obj, method = "ORA", database = "KEGG")
#'
#' # extract enrichment results of downregulated targets
#' enr_down <- targets_enrichment[["downregulated"]]
#'
#' # extract enrichment results as a data.frame
#' enr_df <- enrichmentResults(enr_down)
#'
#' # create a dotplot of enriched terms
#' enrichmentDotplot(enr_down)
#'
#' @note
#' To download gene sets from the above mentioned databases, MIRit uses the
#' `geneset` R package. Moreover, to perform ORA, MIRit implements the
#' `fgsea` package in Bioconductor.
#'
#' @references
#' Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr R
#' package and web server. BMC Bioinformatics 24, 214 (2023).
#' \url{https://doi.org/10.1186/s12859-023-05342-9}.
#' 
#' Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment
#' analysis.” bioRxiv. doi:10.1101/060012,
#' \url{http://biorxiv.org/content/early/2016/06/20/060012}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichTargets <- function(mirnaObj,
                          database = "GO",
                          category = NULL,
                          organism = "Homo sapiens",
                          pCutoff = 0.05,
                          pAdjustment = "fdr",
                          minSize = 10,
                          maxSize = 500) {

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
  if (max(dim(integration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'mirnaIntegration()'",
               "function. See '?mirnaIntegration' for the details."),
         call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("GO", "KEGG", "MsigDB", "WikiPathways", "Reactome",
                       "Enrichr", "DO", "NCG", "DisGeNET", "COVID19")) {
    stop(paste("'database' must be one of 'GO', 'KEGG', 'MsigDB',",
               "'WikiPathways', 'Reactome', 'Enrichr', 'DO', 'NCG',",
               "'DisGeNET', 'COVID19'. For additional details,",
               "see ?enrichTargets"),
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
  if (!is.numeric(minSize) |
      length(minSize) != 1 |
      minSize < 0 |
      !is.numeric(maxSize) |
      length(maxSize) != 1 |
      maxSize < 0 |
      minSize > maxSize) {
    stop(paste("'minSize' and 'maxSize' must be two positive numbers!",
               "Additionally, maxSize should be bigger than minSize.",
               "For additional details, see ?enrichTargets"),
         call. = FALSE)
  }

  ## check if database is supported for the given specie
  supp <- species[!is.na(species[, database]), "specie"]
  if (!organism %in% supp) {
    stop(paste("For", database, "database, 'organism' must be one of:",
               paste(supp, collapse = ", ")))
  }
  
  ## if NULL, set default categories
  if (is.null(category)) {
    if (!database %in% c("Reactome", "WikiPathways", "DO",
                         "DisGeNET", "COVID19")) {
      if (database == "GO") {
        category <- "bp"
      } else if (database == "KEGG") {
        category <- "pathway"
      } else if (database == "MsigDB") {
        category <- "H"
      } else if (database == "NCG") {
        category <- "v7"
      }
      message(paste("Since not specified, 'category' for", database,
                    "database is set to", category, "(default)."))
    } else {
      category <- ""
    }
  }
  
  ## check if category is included in the specified database
  validateCategories(database, category, organism)
  
  ## download the appropriate gene set
  message("Preparing the appropriate gene set...")
  gs <- prepareGeneSet(organism = organism,
                       database = database,
                       category = category)
  
  ## set a describer for the database used
  if (database %in% c("Reactome", "WikiPathways", "DO",
                      "DisGeNET", "COVID19")) {
    dataInfo <- database
  } else {
    dataInfo <- paste(database, " (category: ", category, ")", sep = "")
  }
  
  ## perform over-representation analysis of integrated targets
  res <- oraInternal(mirnaObj = mirnaObj,
                     geneSet = gs,
                     pCutoff = pCutoff,
                     pAdjustment = pAdjustment,
                     minSize = minSize,
                     maxSize = maxSize,
                     dataInfo = dataInfo,
                     organism = organism,
                     integrated = TRUE)
  
  ## return the results of functional enrichment
  return(res)
  
}

