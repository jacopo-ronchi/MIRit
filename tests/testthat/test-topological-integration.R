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
  expect_snapshot(tpRes)
  
})

