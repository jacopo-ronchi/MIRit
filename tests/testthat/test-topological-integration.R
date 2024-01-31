test_that("miRNA augmented pathways are created properly", {
    ## set random seed
    set.seed(1234)

    ## load example objects
    obj <- loadExamples()

    ## prepare augmented pathways
    expect_warning(
        computedPaths <- testPreparePathways(obj),
        "pathways have been ignored because they"
    )

    ## check that every pathway is of class 'graph'
    gEl <- vapply(computedPaths, function(x) is(x, "graph"),
        FUN.VALUE = logical(1)
    )
    allGraphs <- all(gEl == TRUE)
    expect_true(allGraphs)

    ## check the topological sorting
    ths <- computedPaths[["Thyroid hormone synthesis"]]
    expect_snapshot(graph::nodes(ths))
})



test_that("the integrative topological miRNA-mRNA analysis works", {
    ## set random seed
    set.seed(1234)

    ## load example objects
    obj <- loadExamples()
    pt <- loadExamples(class = "IntegrativePathwayAnalysis")

    ## extract the augmented pathways
    paths <- augmentedPathways(pt)

    ## perform topological analysis with TAIPA
    expect_no_error(
        tp <- topologicalAnalysis(obj, paths, nPerm = 100)
    )

    ## extract TAIPA results
    tpRes <- integratedPathways(tp)

    ## check the validity and reproducibility of results
    expect_equal(sum(tpRes$P.Val), 0.019801980198)
})
