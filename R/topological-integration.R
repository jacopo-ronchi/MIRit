topologicalAnalysis <- function(mirnaObj,
                                database = "KEGG",
                                organism = "Homo sapiens",
                                pCutoff = 0.05,
                                pAdjustment = "fdr",
                                nPerm = 10000,
                                minPc = 10,
                                BPPARAM = BiocParallel::bpparam()) {
  
  
  
  ## check inputs
  
  
  
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
  
  ## set organism name and database
  datalow <- tolower(database)
  organism <- convertOrganism(paste("graph", datalow, sep = "_"), organism)
  ## UPDATE !!!!!
  
  ## download pathways from specified database
  message(paste("Downloading pathways from", database, "database..."))
  pathDb <- graphite::pathways(species = organism, database = datalow)
  
  ## convert pathway identifiers to gene symbols
  message("Converting identifiers to gene symbols...")
  pathDb <- quiet(BiocParallel::bplapply(pathDb,
                                         graphite::convertIdentifiers,
                                         to = "SYMBOL",
                                         BPPARAM = BPPARAM))
  ## PREVENT THIS FROM PRINTING !!!!
  
  ## create a list of pathway networks with miRNA-gene interactions
  message("Adding miRNA-gene interactions to biological pathways...")
  netList <- BiocParallel::bplapply(pathDb,
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
  graphList <- BiocParallel::bplapply(netList,
                                      normalizeGraph,
                                      featDE = featDE,
                                      minPc = minPc,
                                      BPPARAM = BPPARAM)
  
  ## remove NULL pathways
  nullPath <- which(lengths(graphList) == 0)
  if (length(nullPath) > 0) {
    warning(paste(length(nullPath), "pathways have been ignored because they",
                  "contain too few nodes with gene expression measurement."))
    graphList <- graphList[-nullPath]
    pNodes <- pNodes[-nullPath]
  }
  
  ## randomly create vector of logFCs with permuted names
  randomVec <- lapply(seq(nPerm), function(i) {
    randomNames <- sample(names(lfc))
    rVec <- setNames(lfc, randomNames)
    rVec
  })
  
  ## compute pathway score and p-value for each pathway
  message(paste("Calculating pathway scores and p-values with",
                nPerm, "permutations..."))
  pS <- lapply_pb(names(graphList), function(pathName) {
    
    ## set graph object and number of nodes
    path <- graphList[[pathName]]
    nN <- pNodes[[pathName]]
    
    ## define pathway score for each graph
    sc <- computePathwayScore(lfc,
                              graph::edges(path),
                              graph::edgeWeights(path))
    
    ## compute pathway score for each random set
    randomScores <- BiocParallel::bplapply(randomVec,
                                           computePathwayScore,
                                           edges = graph::edges(path),
                                           weights = graph::edgeWeights(path),
                                           BPPARAM = BPPARAM)
    randomScores <- unlist(randomScores)
    
    ## define p-value as the fraction of more extreme values in random sets
    pval <- sum(randomScores >= sc) / nPerm
    
    ## calculate normalized score ??? FORSE NON HA SENSO QUI
    norm.sc <- sc / nN
    
    ## return pathway score and p-value for each graph
    c(pathName, sc, norm.sc, pval)
    
  })
  
  ## create result data.frame
  resDf <- do.call(rbind.data.frame, pS)
  colnames(resDf) <- c("pathway", "score", "normalized.score", "P.Val")
  
  ## correct p-values for multiple testing
  resDf$adj.P.Val <- stats::p.adjust(resDf$P.Val, method = pAdjustment)
  
  ## ADD THE NUMBER OF NODES AND THE NUMBER OF MEASURED NODES!!!!
  
  ## order results by adjusted p-values
  resDf <- resDf[order(resDf$adj.P.Val), ]
  
  ## create a TopologicalIntegration object
  res <- new("TopologicalIntegration",
             data = resDf,
             method = paste("Topologically-aware miRNA-mRNA",
                            "integrative analysis (TAMMIA)"),
             organism = organism,
             database = database,
             pCutoff = pCutoff,
             pAdjustment = pAdjustment,
             pathways = netList, ## CORRECT THIS TO INCLUDE NON INVERTED NETS
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
  
  ## determine the number of nodes in the graph
  n <- graph::numNodes(net)
  
  ## remove nodes and edges without measurement
  net <- graph::subGraph(intersect(featDE$ID,
                                   unique(graph::nodes(net))),
                         net)
  
  ## remove nodes that are not connected with any other nodes
  degree <- graph::degree(net)
  dgNode <- degree$inDegree + degree$outDegree
  net <- graph::subGraph(names(dgNode[dgNode != 0]), net)
  
  ## return NULL if graph has less than X % of nodes
  if (graph::numNodes(net) < (minPc / 100) * n) {
    return(NULL)
  }
  
  ## return normalized graph
  return(net)
  
}

