test_that("basic correlation works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform integrative analysis through correlation analysis
  corr <- integrateMirnaTargets(obj)
  
  ## test correlation results
  expect_snapshot(mirnaTargetsIntegration(corr))
  
  ## perform Pearson's correlation analysis with non-default arguments
  corrParam <- integrateMirnaTargets(obj,
                                     pCutoff = 0.03,
                                     pAdjustment = "none",
                                     corMethod = "pearson",
                                     corCutoff = 0.4,
                                     corDirection = "less")
  
  ## test correlation results with non-default arguments
  expect_snapshot(mirnaTargetsIntegration(corrParam))
  
})



test_that("basic Fisher's association works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform integrative analysis through Fisher's exact test
  fsh <- integrateMirnaTargets(obj, test = "fisher")
  
  ## test Fisher's results
  expect_snapshot(mirnaTargetsIntegration(fsh))
  
  ## perform Fisher's association with non-default parameters
  fshParam <- integrateMirnaTargets(obj,
                                    test = "fisher",
                                    pCutoff = 0.2,
                                    pAdjustment = "none",
                                    onlySignificant = FALSE)
  
  ## test Fisher's results with non-default arguments
  expect_snapshot(mirnaTargetsIntegration(fshParam))
  
})
