test_that("NCBI GWAS catalog is reachable through gwasrapidd", {
    ## try to reach the NCBI GWAS catalog
    status <- gwasrapidd::is_ebi_reachable()

    ## check if status is TRUE
    expect_true(status)
})


test_that("Ensembl is reachable through biomaRt", {
    ## skip on windows
    skip_on_os("windows")

    ## try to reach Ensembl
    ensembl <- biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        mirror = "useast"
    )

    ## check if the object is correct
    expect_identical(ensembl@biomart, "ENSEMBL_MART_ENSEMBL")
})
