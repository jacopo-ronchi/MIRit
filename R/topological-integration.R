#' Prepare miRNA-augmented pathways for integrative miRNA-mRNA pathway analyses
#'
#' This function takes influential miRNA-mRNA interactions, identified by the
#' [mirnaIntegration()] function, and adds them to biological pathways
#' retrieved from a pathway database such as `KEGG`, `WikiPathways` and
#' `Reactome`. The pathways returned from this function are needed to perform
#' a topologically-aware integrative pathway analysis (TAIPA) through the
#' [topologicalAnalysis()] function.
#'
#' @details
#' To create augmented pathways, this function uses the `graphite` R package
#' to download biological networks from the above mentioned databases. Then,
#' each pathway is converted to a `graph` object, and significant miRNA-mRNA
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
#' Finally, this function performs a breadth-first search (BFS) algorithm to
#' topologically sort pathway nodes so that each individual node occurs after
#' all its upstream nodes. Nodes within cycles are considered leaf nodes.
#'
#' Information about pathway coverage, i.e. the percentage of nodes with
#' expression measurments, edge weights, topological sorting order, and the
#' parameters used to create the networks are all stored in the `graphData`
#' slot of each `graphNEL` object.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param database The name of the database to use. It must be one of: `KEGG`,
#' `Reactome`, and `WikiPathways`. Default is `KEGG`
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default specie is `Homo sapiens`
#' @param minPc The minimum percentage of measured features that a pathway must
#' have for being considered in the analysis. Default is 10. See the *details*
#' section for additional information
#' @param BPPARAM The desired parallel computing behavior. This parameter
#' defaults to `BiocParallel::bpparam()`, but this can be edited. See
#' [BiocParallel::bpparam()] for information on parallel computing in R
#'
#' @returns
#' A `list` object containing the miRNA-augmented pathways as `graphNEL`
#' objects.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # perform integration analysis with default settings
#' obj <- mirnaIntegration(obj)
#'
#' \donttest{
#' # retrieve pathways from KEGG and augment them with miRNA-gene interactions
#' paths <- preparePathways(obj)
#'
#' # perform the integrative pathway analysis with 1000 permutations
#' ipa <- topologicalAnalysis(obj, paths, nPerm = 1000)
#'
#' # access the results of pathway analysis
#' integratedPathways(ipa)
#'
#' # create a dotplot of integrated pathways
#' integrationDotplot(ipa)
#'
#' # explore a specific biological network
#' visualizeNetwork(ipa, "Thyroid hormone synthesis")
#' }
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
#' @export
preparePathways <- function(mirnaObj,
    database = "KEGG",
    organism = "Homo sapiens",
    minPc = 10,
    BPPARAM = bpparam()) {
    ## input checks
    if (!is(mirnaObj, "MirnaExperiment")) {
        stop("'mirnaObj' should be of class MirnaExperiment! ",
            "See ?MirnaExperiment",
            call. = FALSE
        )
    }
    if (max(dim(integration(mirnaObj))) == 0) {
        stop("Integration analysis is not detected in 'mirnaObj'! ",
            "Before using this function, expression levels of miRNAs and ",
            "genes must be integrated with the 'mirnaIntegration()' ",
            "function. See '?mirnaIntegration' for details.",
            call. = FALSE
        )
    }
    if (!is.character(database) |
        length(database) != 1 |
        !database %in% c("KEGG", "Reactome", "WikiPathways")) {
        stop("Supported databases are: 'KEGG', 'Reactome' and 'WikiPathways'",
            call. = FALSE
        )
    }
    if (!is.numeric(minPc) |
        length(minPc) != 1 |
        minPc < 0 |
        minPc > 100) {
        stop("'minPc' must be a positive number between 0 and 100! ",
            "(default is 10)",
            call. = FALSE
        )
    }

    ## check if database is supported for the given specie
    supp <- species[
        !is.na(species[, paste("graph", database, sep = "_")]),
        "specie"
    ]
    if (!organism %in% supp) {
        stop(
            "For ", database, " database, 'organism' must be one of: ",
            paste(supp, collapse = ", ")
        )
    }

    ## convert organism name to graphite name
    org <- species[
        species$specie == organism,
        paste("graph", database, sep = "_")
    ]

    ## get differential expression results
    deg <- geneDE(mirnaObj, onlySignificant = FALSE)
    dem <- mirnaDE(mirnaObj, onlySignificant = TRUE)

    ## define genes and miRNAs with measurments
    features <- c(dem$ID, deg$ID)

    ## extract integrated miRNA targets
    intTargets <- selectTargets(mirnaObj, pairs = TRUE)

    ## prepare miRNA augmented pathways
    graphList <- preparePathways.internal(
        database = database,
        org = org,
        targ = intTargets,
        features = features,
        minPc = minPc,
        BPPARAM = BPPARAM
    )

    ## prepare graph objects for pathway score calculation
    message("Performing topological sorting of pathway nodes...")
    graphList <- BiocParallel::bplapply(graphList,
        setUpPathways,
        BPPARAM = BPPARAM
    )

    ## add analysis parameters to each graph object
    graphList <- lapply(graphList, function(p) {
        p@graphData$organism <- organism
        p@graphData$database <- database
        p@graphData$minPc <- minPc
        p
    })

    ## returned miRNA-augmented pathways
    return(graphList)
}





#' Perform a topologically-aware integrative pathway analysis (TAIPA)
#'
#' This function allows to perform an integrative pathway analysis that aims
#' to identify the biological networks that are most affected by miRNomic
#' and transcriptomic dysregulations. This function takes miRNA-augmented
#' pathways, created by the [preparePathways()] function, and then calculates
#' a score that estimates the degree of impairment for each pathway. Later,
#' statistical significance is calculated through a permutation test. The main
#' advantages of this method are that it doesn't require matched samples, and
#' that it allows to perform an integrative miRNA-mRNA pathway analysis that
#' takes into account the topology of biological networks. See the *details*
#' section for additional information.
#'
#' @details
#'
#' ## Topologically-Aware Integrative Pathway Analysis (TAIPA)
#'
#' This analysis aims to identify the biological pathways that result affected
#' by miRNA and mRNA dysregulations. In this analysis, biological pathways are
#' retrieved from a pathway database such as KEGG, and the interplay between
#' miRNAs and genes is then added to the networks. Each network is defined as
#' a graph \eqn{G(V, E)}, where \eqn{V} represents nodes, and \eqn{E}
#' represents the relationships between nodes.
#'
#' Then, nodes that are not significantly differentially expressed are assigned
#' a weight \eqn{w_i = 1}, whereas differentially expressed nodes are assigned
#' a weight \eqn{w_i = \left| \Delta E_i \right|}, where \eqn{\Delta E_i} is
#' the linear fold change of the node. Moreover, to consider the biological
#' interaction between two nodes, namely \eqn{i} and \eqn{j}, we define an
#' interaction parameter \eqn{\beta_{i \rightarrow j} = 1} for activation
#' interactions and \eqn{\beta_{i \rightarrow j} = -1} for repression
#' interactions. Subsequently, the concordance coefficient
#' \eqn{\gamma_{i \rightarrow j}} is defined as:
#'
#' \deqn{\gamma_{i \rightarrow j} = \begin{cases} \beta_{i \rightarrow j}
#' &\text{if } sign(\Delta E_i) = sign(\Delta E_j) \\ - \beta_{i \rightarrow j}
#' &\text{if } sign(\Delta E_i) \not= sign(\Delta E_j) \end{cases}\,.}
#'
#' Later in the process, a breadth-first search (BFS) algorithm is applied to
#' topologically sort pathway nodes so that each individual node occurs after
#' all its upstream nodes. Nodes within cycles are considered leaf nodes. At
#' this point, a node score \eqn{\phi} is calculated for each pathway node
#' \eqn{i} as:
#'
#' \deqn{\phi_i = w_i + \sum_{j=1}^{U} \gamma_{i \rightarrow j} \cdot k_j\,.}
#'
#' where \eqn{U} represents the number of upstream nodes,
#' \eqn{\gamma_{i \rightarrow j}} denotes the concordance coefficient, and
#' \eqn{k_j} is a propagation factor defined as:
#'
#' \deqn{k_j = \begin{cases} w_j &\text{if } \phi_j = 0 \\ \phi_j &\text{if }
#' \phi_j \not = 0 \end{cases}\,.}
#'
#' Finally, the pathway score \eqn{\Psi} is calculated as:
#'
#' \deqn{\Psi = \frac{1 - M}{N} \cdot \sum_{i=1}^{N} \phi_i\,,}
#'
#' where \eqn{M} represents the proportion of miRNAs in the pathway, and
#' \eqn{N} represents the total number of nodes in the network.
#'
#' Then, to compute the statistical significance of each pathway score, a
#' permutation procedure is applied. Later, both observed pathway scores and
#' permuted scores are standardized by subtracting the mean score of the
#' permuted sets \eqn{\mu_{\Psi_P}} and then dividing by the standard deviation
#' of the permuted scores \eqn{\sigma_{\Psi_P}}.
#'
#' Finally, the p-value is defined based on the fraction of permutations that
#' reported a higher normalized pathway score than the observed one.
#' However, to prevent p-values equal to zero, we define p-values as:
#'
#' \deqn{p = \frac{\sum_{n=1}^{N_p} \left[ \Psi_{P_N} \ge \Psi_N \right] + 1}
#' {N_p + 1}\,.}
#'
#' In the end, p-values are corrected for multiple testing either through the
#' max-T procedure (default option) which is particularly suited for
#' permutation tests, or through the standard multiple testing approaches.
#'
#' ## Implementation details
#'
#' For computational efficiency, pathway score computation has been implemented
#' in C++ language. Moreover, to define the statistical significance of each
#' network, a permutation test is applied following the number of permutations
#' specified with `nPerm`. The default setting is to perform 10000 permutations.
#' The higher is the number of permutations, the more stable are the calculated
#' p-values, even though the time needed will increase. In this regard, since
#' computing pathway score for 10000 networks for each pathway is
#' computationally intensive, parallel computing has been employed to reduce
#' running time. The user can modify the parallel computing behavior by
#' specifying the `BPPARAM` parameter. See [BiocParallel::bpparam()] for
#' further details. Further, a progress bar can also be included to show the
#' completion percentage by setting `progress = TRUE`. Moreover, the user can
#' define how frequently the progress bar gets updated by tweaking the `tasks`
#' parameter. When using `progress = TRUE`, setting `tasks` to 100 tells the
#' function to update the progress bar 100 times, so that the user can see
#' increases of 1%. Instead, setting `tasks` to 50, means that the progress bar
#' gets updated every 2% of completion. However, keep in mind that `tasks`
#' values from 50 to 100 lead to 15-30% slower p-value calculation due to
#' increased data transfer to the workers. Instead, lower `tasks` values like
#' 20 determine less frequent progress updates but are only slightly less
#' efficient than not including a progress bar.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param pathways A `list` of miRNA-augmented pathways returned by the
#' [preparePathways()] function
#' @param pCutoff The adjusted p-value cutoff to use for statistical
#' significance. The default value is `0.05`
#' @param pAdjustment The p-value correction method for multiple testing. It
#' must be one of: `max-T` (default), `fdr`, `BH`, `none`, `holm`, `hochberg`,
#' `hommel`, `bonferroni`, `BY`
#' @param nPerm The number of permutation used for assessing the statistical
#' significance of each pathway. Default is 10000. See the *details* section
#' for additional information
#' @param progress Logical, whether to show a progress bar during p-value
#' calculation or not. Default is FALSE, not to include a progress bar. Please
#' note that setting `progress = TRUE` with high values of `tasks` leads to
#' less efficient parallelization. See the *details* section for additional
#' information
#' @param tasks An integer between 0 and 100 that specifies how frequently the
#' progress bar must be updated. Default is 0 to simply split the computation
#' among the workers. High values of `tasks` can lead to 15-30% slower p-value
#' calculation. See the *details* section for additional information
#' @param BPPARAM The desired parallel computing behavior. This parameter
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
#' \donttest{
#' # retrieve pathways from KEGG and augment them with miRNA-gene interactions
#' paths <- preparePathways(obj)
#'
#' # perform the integrative pathway analysis with 1000 permutations
#' ipa <- topologicalAnalysis(obj, paths, nPerm = 1000)
#'
#' # access the results of pathway analysis
#' integratedPathways(ipa)
#'
#' # create a dotplot of integrated pathways
#' integrationDotplot(ipa)
#'
#' # explore a specific biological network
#' visualizeNetwork(ipa, "Thyroid hormone synthesis")
#' }
#'
#' @references
#' Peter H. Westfall and S. Stanley Young. Resampling-Based Multiple Testing:
#' Examples and Methods for p-Value Adjustment. John Wiley & Sons.
#' ISBN 978-0-471-55761-6.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
topologicalAnalysis <- function(mirnaObj,
    pathways,
    pCutoff = 0.05,
    pAdjustment = "max-T",
    nPerm = 10000,
    progress = FALSE,
    tasks = 0,
    BPPARAM = bpparam()) {
    ## input checks
    if (!is(mirnaObj, "MirnaExperiment")) {
        stop("'mirnaObj' should be of class MirnaExperiment! ",
            "See ?MirnaExperiment",
            call. = FALSE
        )
    }
    if (max(dim(integration(mirnaObj))) == 0) {
        stop("Integration analysis is not detected in 'mirnaObj'! ",
            "Before using this function, expression levels of miRNAs and ",
            "genes must be integrated with the 'mirnaIntegration()' ",
            "function. See '?mirnaIntegration' for details.",
            call. = FALSE
        )
    }
    if (!is.list(pathways) |
        length(pathways) == 0 |
        !all(vapply(pathways, function(x) is(x, "graph"),
            FUN.VALUE = logical(1)
        ))) {
        stop("'pathways' must be a list object with miRNA-augmented ",
            "pathways returned from the 'preparePathways()' function. ",
            "Please use the 'preparePathways()' function before this call.",
            call. = FALSE
        )
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
            "none", "max-T", "fdr", "bonferroni", "BY",
            "hochberg", "holm", "hommel", "BH"
        )) {
        stop("'pAdjustment' must be  one of: 'none', 'max-T' (default), ",
            "'fdr', 'BH' (same as 'fdr'), 'bonferroni', 'BY', 'hochberg', ",
            "'holm', 'hommel'",
            call. = FALSE
        )
    }
    if (!is.numeric(nPerm) |
        length(nPerm) != 1 |
        nPerm < 0 |
        !nPerm %% 1 == 0) {
        stop("'nPerm' must be a positive integer! (default is 10000)",
            call. = FALSE
        )
    }
    if (!is.logical(progress) |
        length(progress) != 1) {
        stop("'progress' must be logical (TRUE/FALSE)! ",
            "See ?topologicalAnalysis",
            call. = FALSE
        )
    }
    if (!is.numeric(tasks) |
        length(tasks) != 1 |
        tasks < 0 |
        tasks > 100 |
        !tasks %% 1 == 0) {
        stop("'tasks' must be an integer between 0 and 100!",
            call. = FALSE
        )
    }

    ## get differential expression results
    deg <- geneDE(mirnaObj, onlySignificant = FALSE)
    dem <- mirnaDE(mirnaObj, onlySignificant = TRUE)

    ## consider both logFCs and p-values in the definition of weights
    deg$weights <- 1
    sigDeg <- which(deg$ID %in% significantGenes(mirnaObj))
    dem$weights <- 2^abs(dem$logFC)
    deg$weights[sigDeg] <- 2^abs(deg$logFC[sigDeg])

    ## merge miRNA and gene DE results
    dem$type <- "miRNA"
    deg$type <- "gene"
    featDE <- rbind(dem, deg)

    ## calculate the observed score for each pathway
    message("Calculating pathway scores...")
    pS <- bplapply(pathways, function(pathway) {
        computePathwayScore(
            expr = featDE,
            bfs = pathway@graphData$topoSort,
            edges = pathway@graphData$interactions,
            weights = pathway@graphData$eW
        )
    }, BPPARAM = BPPARAM)
    pS <- unlist(pS)

    ## generate n random permutations
    message("Generating random permutations...")
    randomVec <- generatePermutations(deg, dem, nPerm)
    paths <- rep(pathways, times = nPerm)
    permVec <- rep(randomVec, each = length(pathways))

    ## add a functional progress bar for parallel computation
    bpOps <- list()
    if (progress == TRUE) {
        if (hasMethod("bpprogressbar<-",
            signature = c(class(BPPARAM), "logical")
        ) &
            hasMethod("bptasks<-",
                signature = c(class(BPPARAM), "integer")
            )) {
            bpOps <- bpoptions(tasks = tasks, progressbar = TRUE)
        } else {
            warning("Unable to show a progress bar for a BiocParallel backend ",
                "of type ", class(BPPARAM)[1],
                call. = FALSE
            )
        }
    }

    ## use parallel mapply to compute permutation scores for each pathway
    message("Calculating p-values with ", nPerm, " permutations...")
    permScores <- bpmapply(function(pathway, permExpr) {
        computePathwayScore(
            expr = permExpr,
            bfs = pathway@graphData$topoSort,
            edges = pathway@graphData$interactions,
            weights = pathway@graphData$eW
        )
    }, paths, permVec, BPPARAM = BPPARAM, BPOPTIONS = bpOps)

    ## split permuted scores according to pathways
    permList <- split(permScores, names(paths))

    ## order permuted scores based on observed scores
    permList <- permList[names(pS)]

    ## define mean and standard deviation of permuted scores for each pathway
    meanPerm <- lapply(permList, mean)
    stdPerm <- lapply(permList, sd)

    ## normalize observed pathway scores
    normPS <- (as.numeric(pS) - as.numeric(meanPerm)) / as.numeric(stdPerm)
    names(normPS) <- names(pS)

    ## normalize permuted scores as well
    normPerm <- mapply(function(permutedVal, m, sdP) {
        (permutedVal - m) / sdP
    }, permList, meanPerm, stdPerm, SIMPLIFY = FALSE)

    ## define p-values as the fraction of more extreme random values
    pval <- lapply(names(permList), function(pa) {
        (sum(normPerm[[pa]] >= normPS[pa]) + 1) / (nPerm + 1)
    })
    names(pval) <- names(permList)
    pval <- unlist(pval)

    ## extract pathway coverage
    pCov <- vapply(pathways, function(pa) {
        pa@graphData$pathway.coverage
    }, FUN.VALUE = numeric(1))
    names(pCov) <- names(pathways)

    ## create result data.frame
    resDf <- data.frame(
        pathway = names(pS),
        coverage = pCov,
        score = pS,
        normalized.score = normPS,
        P.Val = pval
    )

    ## correct p-values for multiple testing
    if (pAdjustment == "max-T") {
        message("Correcting p-values through the max-T procedure...")
        resDf$adj.P.Val <- maxT.correction(normPerm, normPS, nPerm)
    } else {
        resDf$adj.P.Val <- p.adjust(resDf$P.Val, method = pAdjustment)
    }

    ## order results by adjusted p-values
    resDf <- resDf[order(resDf$adj.P.Val), ]

    ## maintain only significant results
    resDf <- resDf[resDf$adj.P.Val <= pCutoff, ]

    ## print the results of the analysis
    message(
        "The topologically-aware integrative pathway analysis ",
        "reported ", nrow(resDf), " significantly altered pathways!"
    )

    ## create an IntegrativePathwayAnalysis object
    res <- new("IntegrativePathwayAnalysis",
        data = resDf,
        method = paste(
            "Topologically-Aware Integrative",
            "Pathway Analysis (TAIPA)"
        ),
        organism = pathways[[1]]@graphData$organism,
        database = pathways[[1]]@graphData$database,
        pCutoff = pCutoff,
        pAdjustment = pAdjustment,
        pathways = pathways,
        expression = featDE,
        minPc = pathways[[1]]@graphData$minPc,
        nPerm = nPerm
    )

    ## return the results of the integrative analysis
    return(res)
}





## helper function for creating miRNA augmented pathways
preparePathways.internal <- function(database, org, targ, features,
    minPc, BPPARAM) {
  ## load cache
  bfc <- .get_cache()
  
  ## download the appropriate pathways or load them from cache
  dbId <- paste(database, org, sep = "_")
  cache <- BiocFileCache::bfcquery(bfc, dbId)
  if (dbId %in% cache$rname) {
    ## load cached pathways
    message("Reading ", database, " pathways from cache...")
    pathDb <- readRDS(BiocFileCache::bfcrpath(bfc, dbId)[1])
  } else {
    ## download pathways from the specified database
    message("Downloading pathways from ", database, " database...")
    pathDb <- graphite::pathways(species = org, database = tolower(database))
    
    ## retrieve the appropriate organism database
    dbName <- graphite:::selectDb(org)
    
    ## check if the user has installed the required database
    suppressMessages(
      if (!requireNamespace(dbName, quietly = TRUE)) {
        stop("The ", dbName, " package is not installed. Install it ",
             "before runnning this function through: ",
             paste("`BiocManager::install(\"", dbName, "\")`.", sep = ""),
             call. = FALSE
        )
      }
    )
    
    ## convert pathway identifiers to gene symbols by accessing OrgDb through
    ## parallel workers
    message("Converting identifiers to gene symbols...")
    pathDb <- BiocParallel::bplapply(pathDb, function(path) {
      db <- getFromNamespace(dbName, dbName)
      db <- suppressPackageStartupMessages(
        AnnotationDbi::loadDb(AnnotationDbi::dbfile(db))
      )
      on.exit(AnnotationDbi::dbFileDisconnect(AnnotationDbi::dbconn(db)))
      suppressMessages(
        convertNodes(path, db = db)
      )
    }, BPPARAM = BPPARAM)
    
    ## save pathways to cache
    savepath <- BiocFileCache::bfcnew(bfc, rname = dbId, ext = ".RDS")
    saveRDS(pathDb, file = savepath)
  }

    ## create a list of augmented pathways
    message("Adding miRNA-gene interactions to biological pathways...")
    realPaths <- bplapply(pathDb,
        augmentPathway,
        targets = targ,
        BPPARAM = BPPARAM
    )

    ## determine the number of nodes for each augmented pathway
    pNodes <- vapply(realPaths, function(g) {
        if (!is.null(g)) {
            graph::numNodes(g)
        } else {
            NA
        }
    }, FUN.VALUE = integer(1))

    ## normalize networks before computing pathway score
    graphList <- bplapply(realPaths,
        normalizeGraph,
        features = features,
        minPc = minPc,
        BPPARAM = BPPARAM
    )

    ## remove NULL pathways
    nullPath <- which(lengths(graphList) == 0)
    if (length(nullPath) > 0) {
        warning(length(nullPath), " pathways have been ignored because they ",
            "contain too few nodes with gene expression measurement.",
            call. = FALSE
        )
        graphList <- graphList[-nullPath]
    }

    ## obtain the number of considered nodes for each pathway
    cNodes <- vapply(graphList, graph::numNodes, FUN.VALUE = numeric(1))

    ## define pathway coverage as the fraction of considered nodes
    pathCov <- cNodes / pNodes[names(cNodes)]

    ## set pathway coverage as a graph attribute of each pathway
    graphList <- lapply(names(graphList), function(pa) {
        net <- graphList[[pa]]
        net@graphData$pathway.coverage <- pathCov[pa]
        net
    })
    names(graphList) <- names(pathCov)

    ## return pathways
    return(graphList)
}





## helper function for adding miRNA-target pairs to network objects
augmentPathway <- function(pathway,
    targets) {
    ## convert pathway to a graph network
    pathGraph <- graphite::pathwayGraph(pathway = pathway)
    graph::nodes(pathGraph) <- gsub("SYMBOL:", "", graph::nodes(pathGraph))

    ## return NULL if graph does not have nodes
    if (graph::numNodes(pathGraph) == 0) {
        return(NULL)
    }

    ## keep only targets involved in the specified pathway
    pathTargs <- targets[targets$Target %in% graph::nodes(pathGraph), ]

    ## add miRNA-target pairs to the biological network
    pathGraph <- graph::addNode(unique(pathTargs$microRNA), pathGraph)
    pathGraph <- graph::addEdge(pathTargs$microRNA, pathTargs$Target,
        pathGraph,
        weights = -1
    )

    ## assign weights based on interaction type
    ind <- which(as.vector(as(pathGraph, "matrix")) != 0)
    nodeFrom <- (ind - 1) %% length(graph::nodes(pathGraph)) + 1
    nodeTo <- (ind - 1) %/% length(graph::nodes(pathGraph)) + 1
    gMat <- cbind(
        graph::nodes(pathGraph)[nodeFrom],
        graph::nodes(pathGraph)[nodeTo]
    )
    edgeTypes <- unlist(graph::edgeData(
        pathGraph,
        gMat[, 1],
        gMat[, 2],
        "edgeType"
    ))
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

    ## return the graph network
    return(pathGraph)
}





## helper function for removing unmeasured and unconnected nodes
normalizeGraph <- function(net, features, minPc) {
    ## check if graph is NULL
    if (is.null(net)) {
        return(NULL)
    }

    ## determine the number of genes in the graph
    n <- length(graph::nodes(net)[!grepl("miR", graph::nodes(net))])

    ## remove nodes and edges without measurement
    net <- graph::subGraph(
        intersect(
            features,
            unique(graph::nodes(net))
        ),
        net
    )

    ## remove nodes that are not connected with any other nodes
    degree <- graph::degree(net)
    dgNode <- degree$inDegree + degree$outDegree
    net <- graph::subGraph(names(dgNode[dgNode != 0]), net)

    ## return NULL if graph has less than X % of genes with measurement
    nN <- length(graph::nodes(net)[!grepl("miR", graph::nodes(net))])
    if (nN < (minPc / 100) * n |
        nN == 0) {
        return(NULL)
    }

    ## return normalized graph
    return(net)
}





## helper function to topologically sort nodes through a modified BFS algorithm
topologicalSorting <- function(pathway) {
    ## extract node names
    nodes <- graph::nodes(pathway)

    ## initialize levels to -1
    nodeLevels <- rep(-1, length(nodes))
    names(nodeLevels) <- nodes

    ## identify nodes with no ingoing edges (level 0)
    ingEdges <- graph::degree(pathway)$inDegree
    noOut <- names(ingEdges)[ingEdges == 0]
    nodeLevels[noOut] <- 0

    ## calculate levels iteratively
    changed <- TRUE
    currentLevel <- 1
    while (changed) {
        changed <- FALSE
        for (node in nodes) {
            if (nodeLevels[node] == -1) {
                neighborNodes <- unlist(graph::inEdges(node, pathway))
                if (all(nodeLevels[neighborNodes] < currentLevel &
                    all(nodeLevels[neighborNodes] != -1))) {
                    nodeLevels[node] <- currentLevel
                    changed <- TRUE
                }
            }
        }
        currentLevel <- currentLevel + 1
    }

    ## set nodes belonging to cycles to level 0
    nodeLevels[nodeLevels == -1] <- 0

    ## order node levels
    nodeLevels <- nodeLevels[order(nodeLevels)]

    ## return the topological order
    return(nodeLevels)
}





## helper function for setting up graphs as required for score calculation
setUpPathways <- function(pathway) {
    ## perform topological sorting
    pathway@graphData$topoSort <- topologicalSorting(pathway)

    ## create a list with upstream genes for each node
    interactions <- graph::inEdges(graph::nodes(pathway), pathway)
    pathway@graphData$interactions <- interactions

    ## extract edge weights
    eW <- lapply(names(interactions), function(node) {
        nodeInt <- interactions[[node]]
        if (length(nodeInt) > 0) {
            wInt <- graph::edgeData(pathway,
                from = nodeInt,
                to = node, attr = "weight"
            )
            wInt <- unlist(wInt)
            names(wInt) <- nodeInt
            wInt
        } else {
            character(0)
        }
    })
    names(eW) <- names(interactions)

    ## store edge weights in the graph object
    pathway@graphData$eW <- eW

    ## return pathway object
    return(pathway)
}





## helper function for randomly permuting miRNA and gene expression
generatePermutations <- function(deg, dem, nPerm) {
    ## generate n random permutations for miRNAs and genes
    randomDfList <- lapply(seq(nPerm), function(i) {
        ## permute miRNAs
        mirPerm <- dem
        mirPerm$ID <- sample(mirPerm$ID)
        rownames(mirPerm) <- mirPerm$ID

        ## permute genes
        genePerm <- deg
        genePerm$ID <- sample(genePerm$ID)
        rownames(genePerm) <- genePerm$ID

        ## combine miRNAs and genes in a single data.frame
        permDf <- rbind(mirPerm, genePerm)
        permDf
    })

    ## return the computed permutations
    return(randomDfList)
}





## helper function for performing max-T correction on permutation values
maxT.correction <- function(normPerm, normPS, nPerm) {
    ## split permuted scores into a list with each permutation step
    npList <- lapply(seq(nPerm), function(perm) {
        paPerm <- lapply(normPerm, function(pa) {
            pa[perm]
        })
        unlist(paPerm)
    })

    ## record the maximum score for each permutation step
    maxT <- vapply(npList, max, FUN.VALUE = numeric(1))

    ## calculate p-values corrected with max-T procedure
    maxT.pval <- vapply(names(normPS), function(pa) {
        sum(maxT >= normPS[pa]) / nPerm
    }, FUN.VALUE = numeric(1))

    ## return max-T corrected p-values
    return(maxT.pval)
}
