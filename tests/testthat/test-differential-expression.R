test_that("basic differential expression works for edgeR", {
    ## load test object
    obj <- createDummyData()

    ## test differential expression analysis with edgeR
    expect_no_error(
        de <- performMirnaDE(obj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "edgeR"
        )
    )
})


test_that("basic differential expression works for DESeq2", {
    ## skip on Bioconductor
    skip_on_bioc()
    
    ## load test object
    obj <- createDummyData()

    ## test differential expression analysis with DESeq2
    expect_warning(
        de <- performMirnaDE(obj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "DESeq2"
        ),
        "some variables in design formula are characters"
    )
})


test_that("basic differential expression works for limma-voom", {
    ## load test object
    obj <- createDummyData()

    ## test differential expression analysis with limma-voom
    expect_no_error(
        de <- performMirnaDE(obj, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "voom"
        )
    )
})


test_that("basic differential expression works for limma", {
    ## load test object
    obj_microarray <- createDummyData(counts = FALSE)

    ## test differential expression analysis with limma
    expect_no_error(
        de <- performMirnaDE(obj_microarray, "disease", "PTC-NTH",
            ~ 0 + disease,
            method = "limma"
        )
    )
})
