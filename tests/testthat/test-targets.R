test_that("targets retrieval through mirDIP works", {
    ## skip on Bioconductor
    skip_on_bioc()
    
    ## load test object
    obj <- createDummyData(de = TRUE)

    ## retrieve miRNA targets through mirDIP
    expect_no_error(
        obj <- getTargets(obj, includeValidated = FALSE)
    )

    ## check the presence of miRNA-target pairs
    expect_gt(nrow(mirnaTargets(obj)), 0)
})


test_that("miRTarBase is responsive", {
    ## skip on Github Actions and on Bioconductor
    skip_on_ci()
    skip_on_bioc()

    ## targets retrieval through miRTarBase is not tested, but we
    ## check that the miRTarBase link is active
    mtUrl <- paste("https://mirtarbase.cuhk.edu.cn/~miRTarBase/",
        "miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx",
        sep = ""
    )
    hd <- httr::HEAD(mtUrl, httr::timeout(4))
    status <- httr::status_code(hd)
    expect_identical(status, as.integer(200))
})
