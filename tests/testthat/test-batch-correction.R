test_that("batch effect correction works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## correct batch effect
  obj <- batchCorrection(obj, "microRNA", batch = "patient")
  
  ##  verify the resulting matrix
  expect_snapshot(obj[["microRNA"]])
  
})


test_that("batch effect can't be used before differential expression", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## correct batch effect
  obj <- batchCorrection(obj, "microRNA", batch = "patient")
  
  ##  expect warning before running differential expression
  expect_warning(
    obj <- performMirnaDE(obj, "disease", "PTC-NTH",
                          ~ 0 + disease, method = "limma"),
    "Batch effect-corrected matrices"
  )
  
})

