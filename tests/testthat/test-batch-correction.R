test_that("batch effect correction works", {
    ## load the example MirnaExperiment object
    obj <- loadTestObject()

    ## correct batch effect
    expect_no_error(
        obj <- batchCorrection(obj, "microRNA", batch = "patient")
    )

    ##  verify the resulting matrix
    expect_equal(sum(obj[["microRNA"]][, 1]), 370.3815170194314)
    expect_equal(sum(obj[["microRNA"]][, 6]), 397.9221554885439)
})


test_that("batch effect can't be used before differential expression", {
    ## load the example MirnaExperiment object
    obj <- loadTestObject()

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
