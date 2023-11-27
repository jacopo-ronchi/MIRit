test_that("the creation of MirnaExperiment objects works", {
    ## test the creation of the object for paired data
    expect_no_error(
        paired <- createDummyData(paired = "all")
    )

    ## test the creation of the object for partially paired data
    expect_no_error(
        partial <- createDummyData(paired = "partial")
    )

    ## test the creation of the object for unpaired data
    expect_no_error(
        unpaired <- createDummyData(paired = "none")
    )
})
