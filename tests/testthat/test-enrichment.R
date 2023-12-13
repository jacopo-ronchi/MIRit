test_that("gene-sets are correctly created for each supported database", {
    ## check that GO gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "GO", "bp")
    )

    ## check that KEGG gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "KEGG", "pathway")
    )

    ## check that MsigDB gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "MsigDB", "H")
    )

    ## check that WikiPathways gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "WikiPathways", NULL)
    )

    ## check that Reactome gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "Reactome", NULL)
    )

    ## check that Enrichr gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet(
            "Homo sapiens", "Enrichr",
            "Genes_Associated_with_NIH_Grants"
        )
    )

    ## check that DO gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "DO", NULL)
    )

    ## check that NCG gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "NCG", "v7")
    )

    ## check that DisGeNET gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "DisGeNET", NULL)
    )

    ## check that COVID19 gene-sets are correctly retrieved
    expect_no_error(
        gs <- loadSet("Homo sapiens", "COVID19", NULL)
    )
})


test_that("basic ORA enrichment works", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()
    
    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with ORA
    expect_no_error(
        enr <- oraInternal(obj, gs, pCutoff = 0.05, pAdjustment = "none",
                           minSize = 10, maxSize = 500,
                           dataInfo = "KEGG (category: pathway)",
                           organism = "Homo sapiens", integrated = TRUE)
    )

    ## check the validity of ORA
    enrTab <- enrichmentResults(enr$downregulated)
    expect_equal(sum(enrTab$pval), 0.6349798379940161)
})


test_that("basic GSEA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadExamples()
    
    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with GSEA
    expect_no_error(
        enr <- gseaInternal(obj, gs, pCutoff = 0.05, pAdjustment = "none",
                            minSize = 10, maxSize = 500,
                            dataInfo = "KEGG (category: pathway)",
                            organism = "Homo sapiens",
                            rankMetric = "signed.pval", eps = 1e-50)
    )

    ## check the validity of GSEA
    enrTab <- enrichmentResults(enr)
    expect_equal(sum(enrTab$pval), 0.162313766403252)
})


test_that("basic CAMERA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadExamples()
    
    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with CAMERA
    expect_no_error(
        enr <- cameraInternal(obj, gs, pCutoff = 0.05,
                              pAdjustment = "none", minSize = 10,
                              maxSize = 500,
                              dataInfo = "KEGG (category: pathway)",
                              organism = "Homo sapiens")
    )

    ## check the validity of CAMERA
    enrTab <- enrichmentResults(enr)
    expect_equal(sum(enrTab$pval), 0.2732096022757788)
})
