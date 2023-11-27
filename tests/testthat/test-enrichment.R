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

    ## perform functional enrichment with ORA
    expect_no_error(
        enr <- enrichGenes(obj, method = "ORA", database = "GO")
    )

    ## check the validity of ORA
    enrTab <- enrichmentResults(enr$downregulated)
    expect_equal(sum(enrTab$pval), 0.00115495394688)
})


test_that("basic GSEA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadExamples()

    ## perform functional enrichment with GSEA
    expect_no_error(
        enr <- enrichGenes(obj, method = "GSEA", database = "KEGG")
    )

    ## check the validity of GSEA
    enrTab <- enrichmentResults(enr)
    expect_equal(sum(enrTab$pval), 0.00089149223387)
})


test_that("basic CAMERA enrichment works", {
    ## set random seed
    set.seed(1234)

    # load the example MirnaExperiment object
    obj <- loadExamples()

    ## perform functional enrichment with CAMERA
    expect_no_error(
        enr <- enrichGenes(obj,
            method = "CAMERA", database = "KEGG",
            pCutoff = 0.5
        )
    )

    ## check the validity of CAMERA
    enrTab <- enrichmentResults(enr)
    expect_equal(sum(enrTab$pval), 0.0288973773151)
})
