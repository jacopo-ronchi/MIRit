test_that("basic correlation analysis works", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## perform integrative analysis through correlation analysis
    expect_no_error(
        corr <- mirnaIntegration(obj, test = "correlation")
    )

    ## test correlation results
    expect_snapshot(integration(corr))

    ## perform Pearson's correlation analysis with non-default arguments
    expect_no_error(
        corrParam <- mirnaIntegration(obj,
            pCutoff = 0.03,
            pAdjustment = "none",
            corMethod = "pearson",
            corCutoff = 0.4
        )
    )

    ## test correlation results with non-default arguments
    expect_snapshot(integration(corrParam))
})



test_that("basic association tests work", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## perform integrative analysis through Fisher's exact test
    expect_no_error(
        fsh <- mirnaIntegration(obj,
            test = "association",
            associationMethod = "fisher",
            pCutoff = 1
        )
    )

    ## test Fisher's results
    expect_snapshot(integration(fsh))

    ## perform integrative analysis through Fisher's exact test with mid-p
    expect_no_error(
        fshMidP <- mirnaIntegration(obj,
            test = "association",
            associationMethod = "fisher-midp",
            pCutoff = 1
        )
    )

    ## test Fisher's results with mid-p correction
    expect_snapshot(integration(fshMidP))

    ## create a small subset of data to try Boschloo's exact test
    sub <- obj
    sub@mirnaDE$significant <- sub@mirnaDE$significant[seq(2)]

    ## test the integration through Boschloo's test for just two miRNAs
    expect_no_error(
        bosch <- mirnaIntegration(sub,
            test = "association",
            associationMethod = "boschloo",
            pCutoff = 1
        )
    )

    ## test Boschloo's results
    expect_snapshot(integration(bosch))
})


test_that("basic fry integration works", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## perform integrative analysis through fry
    expect_no_error(
        fry <- mirnaIntegration(obj, test = "fry", pCutoff = 1)
    )

    ## test fry results
    expect_snapshot(integration(fry))
})


test_that("correlation analysis works for partially paired datasets", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## change object metadata to define unpaired samples
    unp <- createDummyData(
        nGenes = nrow(obj[["genes"]]),
        nMirnas = nrow(obj[["microRNA"]]),
        counts = FALSE, paired = "partial"
    )
    unp@mirnaDE <- obj@mirnaDE
    unp@geneDE <- obj@geneDE
    unp@targets <- obj@targets

    ## perform correlation analysis on partially matched samples
    expect_warning(
        corr <- mirnaIntegration(unp, test = "correlation"),
        "Removed samples are: PTC 7, PTC 8, NTH 7, NTH 8, ptc_7, ptc_8,"
    )

    ## test correlation results
    expect_snapshot(integration(corr))
})
