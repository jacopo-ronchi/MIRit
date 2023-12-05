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





## helper function to test the retrieval of gene-sets
loadSet <- function(organism, database, category) {
    ## determine organism name accepted by database
    org <- species[species$specie == organism, database]

    ## determine organism name for ID conversion
    convOrg <- species[species$specie == organism, "Conversion"]

    ## download and prepare the appropriate gene set
    if (database == "GO") {
        gs <- geneset::getGO(
            org = org,
            ont = category
        )
    } else if (database == "KEGG") {
        gs <- geneset::getKEGG(
            org = org,
            category = category
        )
    } else if (database == "MsigDB") {
        gs <- geneset::getMsigdb(
            org = org,
            category = category
        )
    } else if (database == "WikiPathways") {
        gs <- geneset::getWiki(org = org)
    } else if (database == "Reactome") {
        gs <- geneset::getReactome(org = org)
    } else if (database == "Enrichr") {
        gs <- geneset::getEnrichrdb(
            org = org,
            library = category
        )
    } else if (database == "DO") {
        gs <- geneset::getHgDisease(source = "do")
    } else if (database == "NCG") {
        gs <- geneset::getHgDisease(source = "ncg_v6")
    } else if (database == "DisGeNET") {
        gs <- geneset::getHgDisease(source = "disgenet")
    } else if (database == "COVID19") {
        gs <- geneset::getHgDisease(source = "covid19")
    }

    ## return gene set
    return(gs)
}
