test_that("the creation of MirnaExperiment objects works", {
    ## test the creation of the object for paired data
    paired <- createDummyData(paired = "all")
    expect_snapshot(paired)

    ## test the creation of the object for partially paired data
    partial <- createDummyData(paired = "partial")
    expect_snapshot(partial)

    ## test the creation of the object for unpaired data
    unpaired <- createDummyData(paired = "none")
    expect_snapshot(unpaired)
})
