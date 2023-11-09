test_that("targets retrieval through mirDIP works", {
  
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
  
  ## targets retrieval through miRTarBase is not tested, but we
  ## check that the miRTarBase link is active
  mtUrl <- paste("https://mirtarbase.cuhk.edu.cn/~miRTarBase/",
                 "miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx",
                 sep = "")
  hd <- httr::HEAD(mtUrl)
  status <- hd$all_headers[[1]]$status
  expect_identical(status, as.integer(200))
  
})

