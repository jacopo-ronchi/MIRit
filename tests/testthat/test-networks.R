test_that("miRNA-mRNA network construction works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## create miRNA-mRNA KEGG network
  pKegg <- mirnaPathway(obj,
                        pathway = "Thyroid hormone synthesis",
                        database = "KEGG",
                        organism = "Homo sapiens")
  
  ## test KEGG network results
  expect_snapshot(pKegg)
  
  ## create miRNA-mRNA Reactome network
  pReac <- mirnaPathway(obj,
                        pathway = "Regulation of thyroid hormone activity",
                        database = "Reactome",
                        organism = "Homo sapiens")
  
  ## test Reactome network results
  expect_snapshot(pReac)
  
})
