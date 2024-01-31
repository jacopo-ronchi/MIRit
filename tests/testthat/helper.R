## create dummy objects needed for running unit tests
createDummyData <- function(
        nGenes = 100,
        nMirnas = 50,
        counts = TRUE,
        paired = "all",
        de = FALSE) {
    ## load experiment
    obj <- loadExamples()

    ## extract miRNA and gene expression measurements
    if (counts == TRUE) {
        geneExpr <- obj@metadata$oldCounts$genes
        mirnaExpr <- obj@metadata$oldCounts$microRNA
    } else {
        geneExpr <- obj[["genes"]]
        mirnaExpr <- obj[["microRNA"]]
    }

    ## resize expression data
    geneExpr <- geneExpr[seq(nGenes), ]
    mirnaExpr <- mirnaExpr[seq(nMirnas), ]

    ## define paired and unpaired samples
    if (paired == "all") {
        pSamp <- TRUE
        meta <- data.frame(
            primary = colnames(mirnaExpr),
            mirnaCol = colnames(mirnaExpr),
            geneCol = colnames(geneExpr),
            disease = c(rep("PTC", 8), rep("NTH", 8))
        )
    } else if (paired == "partial") {
        pSamp <- TRUE
        extraSamples <- c("ptc_7", "ptc_8", "nth_7", "nth_8")
        colnames(geneExpr)[c(7, 8, 15, 16)] <- extraSamples
        meta <- data.frame(
            primary = c(colnames(mirnaExpr), extraSamples),
            mirnaCol = c(colnames(mirnaExpr), rep(NA, 4)),
            geneCol = c(
                colnames(geneExpr)[seq(6)], NA, NA,
                colnames(geneExpr)[seq(9, 14)], NA, NA,
                extraSamples
            ),
            disease = c(
                rep("PTC", 8), rep("NTH", 8),
                rep("PTC", 2), rep("NTH", 2)
            )
        )
    } else if (paired == "none") {
        pSamp <- FALSE
        colnames(geneExpr) <- c(
            paste("ptc", seq(8), sep = "_"),
            paste("nth", seq(8), sep = "_")
        )
        meta <- data.frame(
            primary = c(colnames(mirnaExpr), colnames(geneExpr)),
            mirnaCol = c(colnames(mirnaExpr), rep(NA, 16)),
            geneCol = c(rep(NA, 16), colnames(geneExpr)),
            disease = c(
                rep("PTC", 8), rep("NTH", 8),
                rep("PTC", 8), rep("NTH", 8)
            )
        )
    }

    ## create a dummy MirnaExperiment object
    dummyObj <- MirnaExperiment(mirnaExpr, geneExpr,
        samplesMetadata = meta,
        pairedSamples = pSamp
    )

    ## perform differential expression analysis
    if (de == TRUE) {
        dummyObj <- performMirnaDE(dummyObj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "edgeR",
            pAdjustment = "none"
        )
        dummyObj <- performGeneDE(dummyObj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "edgeR",
            pAdjustment = "none"
        )
    }

    ## return dummy object
    return(dummyObj)
}





## helper function to load KEGG gene-sets for testing enrichment methods
loadExampleGeneSets <- function() {
    ## load the example GSEA object
    enr <- loadExamples("FunctionalEnrichment")
    
    ## extract gene-sets from example GSEA object
    gs <- geneSet(enr)
    
    ## return the gene-sets
    return(gs)
}





## helper function for testing the construction of augmented pathways
testPreparePathways <- function(
        mirnaObj,
        database = "KEGG",
        organism = "Homo sapiens",
        minPc = 10,
        size = NULL,
        BPPARAM = bpparam()) {
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
    graphList <- testPreparePathways.internal(
        database = database,
        org = org,
        targ = intTargets,
        features = features,
        minPc = minPc,
        size = size,
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
    
    ## return miRNA-augmented pathways
    return(graphList)
}





## helper function to prepare augmented pathways from a subset of KEGG pathways
testPreparePathways.internal <- function(
        database, org, targ, features,
        minPc, size, BPPARAM) {
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
        pathDb <- graphite::pathways(
            species = org,
            database = tolower(database)
        )
        
        ## retrieve the appropriate organism database
        dbName <- graphite:::selectDb(org)
        
        ## check if the user has installed the required database
        suppressMessages(
            if (!requireNamespace(dbName, quietly = TRUE)) {
                stop("The ", dbName, " package is not installed. Install it ",
                     "before runnning this function through: ",
                     paste("`BiocManager::install(\"", dbName, "\")`.",
                           sep = ""
                     ),
                     call. = FALSE
                )
            }
        )
        
        ## convert pathway identifiers to gene symbols by accessing OrgDb
        ## through parallel workers
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
    
    ## select a subset of KEGG pathways for testing purposes
    pathDb <- pathDb[205:210]
    
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
    
    ## remove pathways with too many nodes or too few nodes
    if (!is.null(size)) {
        graphList <- graphList[cNodes >= size[1] & cNodes <= size[2]]
    }
    
    ## return pathways
    return(graphList)
}
