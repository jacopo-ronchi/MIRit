test_that("batch effect correction works", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## correct batch effect
    expect_no_error(
        obj <- batchCorrection(obj, "microRNA", batch = "patient")
    )

    ##  verify the resulting matrix
    expect_equal(sum(obj[["microRNA"]][, 1]), 2529.95084387)
    expect_equal(sum(obj[["microRNA"]][, 6]), 2542.02835449)
})


test_that("batch effect can't be used before differential expression", {
    ## load the example MirnaExperiment object
    obj <- loadExamples()

    ## correct batch effect
    obj <- batchCorrection(obj, "microRNA", batch = "patient")

    ##  expect warning before running differential expression
    expect_warning(
        obj <- performMirnaDE(obj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "limma"
        ),
        "Batch effect-corrected matrices"
    )
})
