test_that("basic differential expression works for edgeR", {
  
  ## load test object
  obj <- createDummyData()
  
  ## test differential expression analysis with edgeR
  de_edgeR <- performMirnaDE(obj, "disease", "PTC-NTH",
                             ~ 0 + disease, method = "edgeR")
  expect_snapshot(mirnaDE(de_edgeR))
  
})


test_that("basic differential expression works for DESeq2", {
  
  ## load test object
  obj <- createDummyData()
  
  ## test differential expression analysis with DESeq2
  de_DESeq2 <- performMirnaDE(obj, "disease", "PTC-NTH",
                              ~ 0 + disease, method = "DESeq2")
  expect_snapshot(mirnaDE(de_DESeq2))
  
})


test_that("basic differential expression works for limma-voom", {
  
  ## load test object
  obj <- createDummyData()
  
  ## test differential expression analysis with limma-voom
  de_voom <- performMirnaDE(obj, "disease", "PTC-NTH",
                            ~ 0 + disease, method = "voom")
  expect_snapshot(mirnaDE(de_voom))
  
})


test_that("basic differential expression works for limma", {
  
  ## load test object
  obj_microarray <- createDummyData(counts = FALSE)
  
  ## test differential expression analysis with limma
  de_limma <- performMirnaDE(obj_microarray, "disease", "PTC-NTH",
                             ~ 0 + disease, method = "limma")
  expect_snapshot(mirnaDE(de_limma))
  
})

