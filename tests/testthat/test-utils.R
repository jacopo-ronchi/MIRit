test_that("integrated targets selection works for each integration strategy", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## integrate data through the various approaches (except for correlation)
    asc <- mirnaIntegration(obj,
        test = "association",
        associationMethod = "fisher",
        pCutoff = 0.5,
        pAdjustment = "none"
    )
    fry <- mirnaIntegration(obj,
        test = "fry",
        pCutoff = 0.5,
        pAdjustment = "none"
    )

    ## test integrated pairs extraction for correlation analysis
    inTr <- selectTargets(obj, pairs = TRUE)
    expect_snapshot(inTr)

    ## test integrated pairs extraction for asscoiation tests
    inTrAsc <- selectTargets(asc, pairs = TRUE)
    expect_snapshot(inTrAsc)

    ## test integrated pairs extraction for fry method
    inTrFry <- selectTargets(fry, pairs = TRUE)
    expect_snapshot(inTrFry)
})
