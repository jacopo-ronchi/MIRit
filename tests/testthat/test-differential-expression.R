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
    deTab <- mirnaDE(de)
    expect_equal(sum(deTab$P.Value), 0.00812718819306)
    expect_equal(sum(deTab$logFC), 1.45879204472)
})


test_that("basic differential expression works for DESeq2", {
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
    deTab <- mirnaDE(de)
    expect_equal(sum(deTab$P.Value), 0.000296432857819)
    expect_equal(sum(deTab$logFC), -2.07588192405)
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
    deTab <- mirnaDE(de)
    expect_equal(sum(deTab$P.Value), 0.000360018745351)
    expect_equal(sum(deTab$logFC), -1.05568453494)
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
    deTab <- mirnaDE(de)
    expect_equal(sum(deTab$P.Value), 0.000918172142831)
    expect_equal(sum(deTab$logFC), -3.50564733182)
})
