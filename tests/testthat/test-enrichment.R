test_that("basic ORA enrichment works", {
    ## load the example MirnaExperiment object
    obj <- loadTestObject()

    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with ORA
    expect_no_error(
        enr <- oraInternal(obj, gs,
            pCutoff = 0.05, pAdjustment = "none",
            minSize = 0, maxSize = 500,
            dataInfo = "KEGG (category: pathway)",
            organism = "Homo sapiens", integrated = TRUE
        )
    )

    ## check the validity of ORA
    enrTab <- enrichmentResults(enr$downregulated)
    expect_equal(sum(enrTab$pval), 0.03193418836636841)
})


test_that("basic GSEA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadTestObject()

    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with GSEA
    expect_no_error(
        enr <- gseaInternal(obj, gs,
            pCutoff = 0.4, pAdjustment = "none",
            minSize = 0, maxSize = 500,
            dataInfo = "KEGG (category: pathway)",
            organism = "Homo sapiens",
            rankMetric = "signed.pval", eps = 1e-50
        )
    )
})


test_that("basic CAMERA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadTestObject()

    ## load example gene-sets
    gs <- loadExampleGeneSets()

    ## perform functional enrichment with CAMERA
    expect_no_error(
        enr <- cameraInternal(obj, gs,
            pCutoff = 0.1,
            pAdjustment = "none", minSize = 10,
            maxSize = 500,
            dataInfo = "KEGG (category: pathway)",
            organism = "Homo sapiens"
        )
    )

    ## check the validity of CAMERA
    enrTab <- enrichmentResults(enr)
    expect_equal(sum(enrTab$pval), 0.06110395732354543)
})
