#' Perform a topologically-aware integrative pathway analysis (TAIPA)
#'
#' This function allows to perform an integrative pathway analysis that aims
#' to identify the biological networks that are most affected by miRNomic
#' and transcriptomic dysregulations. Briefly, influential miRNA-mRNA
#' interactions, identified by the [mirnaIntegration()] function, are added to
#' biological pathways retrieved from a pathway database such as `KEGG`,
#' `WikiPathways` and `Reactome`. Then, a score that estimates the degree of
#' impairment is calculated for each pathway, and statistical significance is
#' calculated through a permutation test. The main advantages of this method
#' is that it doesn't require matched samples, and that it allows to perform
#' an integrative miRNA-mRNA pathway analysis that take into account the
#' topology of biological pathways and miRNA-mRNA interactions. See the
#' *details* section for additional information.
#'
#' @details
#' 
#' ## Topologically-Aware Integrative Pathway Analysis (TAIPA)
#' 
#' This analysis aims to identify the biological pathways that result affected
#' by miRNA and mRNA dysregulations. In this analysis, biological pathways are
#' retrieved from a pathway database such as KEGG, and the interplay between
#' miRNAs and genes is then added to the networks. Each network is defined as
#' a graph \eqn{G(N,E)}, where \eqn{N} represents nodes, and \eqn{E} represents
#' the relationships between nodes. For each pathway, a score is a computed as:
#' 
#' \deqn{s = \frac{1}{N}\cdot \sum_{i = 0}^{N}\sum_{j = 0}^{U}\beta_{i\to j}
#' \left( \Delta E_i\cdot \Delta E_j\right)\,,}
#' 
#' where \eqn{N} is the total number of nodes in the pathway, \eqn{U} is the
#' total number of interactors for each node, \eqn{\Delta E} is the log2 fold
#' change of each node, and \eqn{\beta} is an interaction factor that
#' represents the direction of a specific interaction (+1 for activation,
#' -1 for inhibition).
#' 
#' The higher is the score, and the more affected is a pathway. Then, to
#' compute the statistical significance of each pathway score, a permutation
#' procedure is applied. Indeed, for each pathway log2 fold changes are
#' permuted several times, and for each permutation pathway score is calculated.
#' In the end, the p-value is defined as the fraction of permutations that
#' reported a higher pathway score than the observed one.
#' 
#' ## Implementation details
#' 
#' To create augmented pathways, this function uses the `graphite` R package
#' to download biological networks from the above mentioned database. Then,
#' each pathway is converted to a graph object, and significant miRNA-mRNA
#' interactions are added to the network. Further, edge weights are added
#' according to interaction type.
#' 
#' At this point, biological pathways with few nodes measured are excluded
#' from this analysis. This is required because, during differential expression
#' analysis, lowly expressed features are removed. Therefore, some pathways
#' might result significantly affected even if only 1% of nodes is perturbed.
#' The default behavior is to exclude pathways with less than 10% of
#' representation (`minPc = 10`).
#' 
#' After this normalization step, the score of each pathway is calculated.
#' For computational efficiency, pathway score computation has been implemented
#' in C++ language.
#' 
#' Moreover, to define the statistical significance of each network, a
#' permutation test is applied following the number of permutations specified
#' with `nPerm`. The default setting is to perform 10000 permutations. The
#' higher is the number of permutations, the more stable are the calculated
#' p-values, even though the time needed will increase. In this regard, since
#' computing pathway score for 10000 networks for each pathway is
#' computationally intensive, parallel computing has been employed to reduce
#' running time. The user can modify the parallel computing behavior by
#' specifying the `BPPARAM` parameter. See [BiocParallel::bpparam()] for
#' further details.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database to use. It must be one of: `KEGG`,
#' `Reactome`, and `WikiPathways`. Default is `KEGG`
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default specie is `Homo sapiens`
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
#' `bonferroni`, `BY`
#' @param nPerm The number of permutation used for assessing the statistical
#' significance of each pathway. Default is 10000. See the *details* section
#' for additional information
#' @param minPc The minimum percentage of measured features that a pathway must
#' have for being considered in the analysis. Default is 10. See the *details*
#' section for additional information
#' @param BPPARAM The desired parallel computing behavior. This paramete
#' defaults to `BiocParallel::bpparam()`, but this can be edited. See
#' [BiocParallel::bpparam()] for information on parallel computing in R
#'
#' @returns
#' An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] that stores
#' the results of the analysis. See the relative help page for further details.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform integration analysis with default settings
#' obj <- mirnaIntegration(obj)
#'
#' # perform the integrative pathway analysis with default settings
#' #ipa <- topologicalAnalysis(obj)
#' 
#' # access the results of pathway analysis
#' #integratedPathways(ipa)
#' 
#' # create a dotplot of integrated pathways
#' #integrationDotplot(ipa)
#' 
#' # explore a specific biological network
#' #visualizeNetwork(ipa, "Thyroid hormone synthesis")
#' 
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012),
#' \url{https://doi.org/10.1186/1471-2105-13-20}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#' 
#' @import BiocParallel
#' @export
topologicalAnalysis <- function(mirnaObj,
                                database = "KEGG",
                                organism = "Homo sapiens",
                                pCutoff = 0.05,
                                pAdjustment = "fdr",
                                nPerm = 10000,
                                minPc = 10,
                                BPPARAM = bpparam()) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (max(dim(integration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'mirnaIntegration()'",
               "function. See '?mirnaIntegration' for details."),
         call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("KEGG", "Reactome", "WikiPathways")) {
    stop("Supported databases are: 'KEGG', 'Reactome' and 'WikiPathways'",
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
  if (!is.numeric(nPerm) |
      length(nPerm) != 1 |
      nPerm < 0 |
      !nPerm %% 1 == 0) {
    stop("'nPerm' must be a positive integer! (default is 10000)",
         call. = FALSE)
  }
  if (!is.numeric(minPc) |
      length(minPc) != 1 |
      minPc < 0 |
      minPc > 100) {
    stop("'minPc' must be a positive number between 0 and 100! (default is 10)",
         call. = FALSE)
  }
  
  ## check if database is supported for the given specie
  supp <- species[!is.na(species[, paste("graph", database, sep = "_")]),
                  "specie"]
  if (!organism %in% supp) {
    stop(paste("For", database, "database, 'organism' must be one of:",
               paste(supp, collapse = ", ")))
  }
  
  ## convert organism name to graphite name
  org <- species[species$specie == organism,
                 paste("graph", database, sep = "_")]
  
  ## get differential expression results
  deg <- geneDE(mirnaObj, onlySignificant = FALSE)
  dem <- mirnaDE(mirnaObj, onlySignificant = FALSE)
  featDE <- rbind(dem, deg)
  
  ## create a named vector with logFCs for each feature
  lfc <- featDE$logFC
  names(lfc) <- featDE$ID
  
  ## extract integrated miRNA targets
  intTargets <- c(selectTargets(mirnaObj, "upregulated"),
                  selectTargets(mirnaObj, "downregulated"))
  targ <- mirnaTargets(mirnaObj)
  targ <- targ[targ$Gene.Symbol %in% intTargets, ]
  
  ## download pathways from specified database
  message(paste("Downloading pathways from", database, "database..."))
  pathDb <- graphite::pathways(species = org, database = tolower(database))
  
  ## convert pathway identifiers to gene symbols
  message("Converting identifiers to gene symbols...")
  pathDb <- bplapply(as.list(pathDb), function(x) {
    suppressMessages(graphite::convertIdentifiers(x, to = "SYMBOL"))
  }, BPPARAM = BPPARAM)
  
  ## create a list of augmented pathways without inversion
  message("Adding miRNA-gene interactions to biological pathways...")
  realPaths <- bplapply(pathDb,
                        augmentPathway,
                        targets = targ,
                        inverted = FALSE,
                        BPPARAM = BPPARAM)
  
  ## create a list of inverted networks for pathway score calculation
  netList <- bplapply(pathDb,
                      augmentPathway,
                      targets = targ,
                      inverted = TRUE,
                      BPPARAM = BPPARAM)
  
  ## determine the number of nodes for each augmented pathway
  pNodes <- vapply(netList, function(g) {
    if (!is.null(g)) {
      graph::numNodes(g)
    } else {
      NA
    }
  }, FUN.VALUE = integer(1))
  
  ## normalize networks before computing pathway score
  graphList <- bplapply(netList,
                        normalizeGraph,
                        featDE = featDE,
                        minPc = minPc,
                        BPPARAM = BPPARAM)
  
  ## remove NULL pathways
  nullPath <- which(lengths(graphList) == 0)
  if (length(nullPath) > 0) {
    warning(paste(length(nullPath), "pathways have been ignored because they",
                  "contain too few nodes with gene expression measurement."),
            call. = FALSE)
    graphList <- graphList[-nullPath]
    pNodes <- pNodes[-nullPath]
  }
  
  ## obtain the number of considered nodes for each pathway
  cNodes <- unlist(lapply(graphList, graph::numNodes))
  
  ## randomly create vector of logFCs with permuted names
  randomVec <- lapply(seq(nPerm), function(i) {
    randomNames <- sample(names(lfc))
    rVec <- setNames(lfc, randomNames)
    rVec
  })
  
  ## calculate the observed score for each pathway
  message("Calculating pathway scores...")
  pS <- bplapply(names(graphList),
                 score.helper,
                 lfc = lfc,
                 gList = graphList,
                 pN = pNodes,
                 BPPARAM = BPPARAM)
  names(pS) <- names(graphList)
  pS <- unlist(pS)
  
  ## combine pathway networks and permuted vectors
  paths <- rep(names(graphList), times = nPerm)
  permVec <- rep(randomVec, each = length(graphList))
  
  ## use parallel mapply to compute score for each pathway and for each set
  message(paste("Calculating p-values with", nPerm, "permutations..."))
  permScores <- bpmapply(score.helper,
                         paths,
                         permVec,
                         MoreArgs = list(gList = graphList,
                                         pN = pNodes),
                         BPPARAM = BPPARAM,
                         BPOPTIONS = bpoptions(tasks = 500,
                                               progressbar = TRUE))
  
  ## split the permuted scores according to the pathways
  permScores <- split(permScores, paths)
  
  ## order permuted scores based on observed scores
  permScores <- permScores[names(pS)]
  
  ## define p-values as the fraction of more extreme random values
  pval <- lapply(names(permScores), function(pa) {
    sum(permScores[[pa]] >= pS[[pa]]) / nPerm
  })
  names(pval) <- names(permScores)
  pval <- unlist(pval)
  
  ## create result data.frame
  resDf <- data.frame(pathway = names(pS),
                      considered.nodes = cNodes,
                      total.nodes = pNodes,
                      score = pS,
                      P.Val = pval)
  
  ## correct p-values for multiple testing
  resDf$adj.P.Val <- stats::p.adjust(resDf$P.Val, method = pAdjustment)
  
  ## order results by adjusted p-values
  resDf <- resDf[order(resDf$adj.P.Val), ]
  
  ## maintain only significant results
  resDf <- resDf[resDf$adj.P.Val <= pCutoff, ]
  
  ## print the results of the analysis
  message(paste("The topologically-aware integrative pathway analysis",
                "reported", nrow(resDf), "significantly altered pathways!"))
  
  ## create a IntegrativePathwayAnalysis object
  res <- new("IntegrativePathwayAnalysis",
             data = resDf,
             method = paste("Topologically-Aware Integrative",
                            "Pathway Analysis (TAIPA)"),
             organism = organism,
             database = database,
             pCutoff = pCutoff,
             pAdjustment = pAdjustment,
             pathways = realPaths,
             expression = lfc,
             minPc = minPc,
             nPerm = nPerm)
  
  ## return the results of the integrative analysis
  return(res)
  
}





## helper function for creating and normalizing miRNA augmented pathways
augmentPathway <- function(pathway,
                           targets,
                           inverted = FALSE) {
  
  ## convert pathway to a graph network
  pathGraph <- graphite::pathwayGraph(pathway = pathway)
  graph::nodes(pathGraph) <- gsub("SYMBOL:", "", graph::nodes(pathGraph))
  
  ## return NULL if graph does not have nodes
  if (graph::numNodes(pathGraph) == 0) {
    return(NULL)
  }
  
  ## keep only targets involved in the specified pathway
  pathTargs <- targets[targets$Gene.Symbol %in% graph::nodes(pathGraph), ]
  
  ## add miRNA-target pairs to the biological network
  pathGraph <- graph::addNode(unique(pathTargs$MicroRNA), pathGraph)
  pathGraph <- graph::addEdge(pathTargs$MicroRNA, pathTargs$Gene.Symbol,
                              pathGraph, weights = -1)
  
  ## assign weights based on interaction type
  ind <- which(as.vector(as(pathGraph, "matrix")) != 0)
  nodeFrom <- (ind - 1) %% length(graph::nodes(pathGraph)) + 1
  nodeTo <- (ind - 1) %/% length(graph::nodes(pathGraph)) + 1
  gMat <- cbind(graph::nodes(pathGraph)[nodeFrom],
                graph::nodes(pathGraph)[nodeTo])
  edgeTypes <- unlist(graph::edgeData(pathGraph,
                                      gMat[, 1],
                                      gMat[, 2],
                                      "edgeType"))
  edgeTypes[edgeTypes == "undefined"] <- "inhibition"
  
  w <- vapply(edgeTypes, function(type) {
    activation <- "activation|expression"
    inhibition <- "inhibition|repression"
    if (grepl(activation, type) & grepl(inhibition, type) == TRUE) {
      assignedWeight <- 0
    } else if (grepl(activation, type) == TRUE) {
      assignedWeight <- 1
    } else if (grepl(inhibition, type) == TRUE) {
      assignedWeight <- -1
    } else {
      assignedWeight <- 0
    }
  }, FUN.VALUE = numeric(1))
  
  suppressWarnings(
    pathGraph <- graph::addEdge(gMat[, 1], gMat[, 2], pathGraph, w)
  )
  
  ## invert graph object for pathway score calculation
  if (inverted == TRUE) {
    pathGraph <- graph::reverseEdgeDirections(pathGraph)
    suppressWarnings(
      pathGraph <- graph::addEdge(gMat[, 2], gMat[, 1], pathGraph, w)
    )
  }
  
  ## return the graph network
  return(pathGraph)
  
}





## helper function for removing unmeasured and unconnected nodes
normalizeGraph <- function(net, featDE, minPc) {
  
  ## check if graph is NULL
  if (is.null(net)) {
    return(NULL)
  }
  
  ## determine the number of genes in the graph
  n <- length(graph::nodes(net)[!grepl("miR", graph::nodes(net))])
  
  ## remove nodes and edges without measurement
  net <- graph::subGraph(intersect(featDE$ID,
                                   unique(graph::nodes(net))),
                         net)
  
  ## remove nodes that are not connected with any other nodes
  degree <- graph::degree(net)
  dgNode <- degree$inDegree + degree$outDegree
  net <- graph::subGraph(names(dgNode[dgNode != 0]), net)
  
  ## return NULL if graph has less than X % of genes with measurement
  nN <- length(graph::nodes(net)[!grepl("miR", graph::nodes(net))])
  if (nN < (minPc / 100) * n) {
    return(NULL)
  }
  
  ## return normalized graph
  return(net)
  
}





## helper function for calculating pathway score for each graph in a list
score.helper <- function(pathName, lfc, gList, pN) {
  
  ## set graph object and number of nodes
  path <- gList[[pathName]]
  nN <- pN[[pathName]]
  
  ## define pathway score for each graph
  sc <- computePathwayScore(lfc,
                            edges = graph::edges(path),
                            weights = graph::edgeWeights(path),
                            nN = nN)
  
  ## return the resulting score
  return(sc)
  
}

