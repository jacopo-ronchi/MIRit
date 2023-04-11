test_that("listCategories works", {
  
  ## retrieve available categories in miEAA 2.0
  available_categories <- listCategories(organism = "Mus musculus",
                                         mirnaType = "precursor")
  
  ## define expected categories
  expected <- c(`Chromosomal location (miRBase)` =
                  "miRBase_Chromosomes_precursor", 
                `Cluster (miRBase)` = "miRBase_Cluster_precursor",
                `Family (miRBase)` = "miRBase_Family_precursor", 
                `Confidence (miRBase)` = "miRBase_High_confidence_precursor", 
                `Pubmed (miRBase)` = "miRBase_Pubmed_precursor",
                `Diseases (MNDR)` = "MNDR_precursor", 
                `Interactions (NPInter)` = "NPInter_precursor",
                `Expressed in tissue (Atlas)` =
                  "Published_expressed_tissue_precursor", 
                `Tissue specific (Atlas)` =
                  "Published_specific_tissue_precursor", 
                `Localization (RNALocate)` = "RNALocate_precursor",
                `Tissues with miRNA-TF interactions (TransmiR)` =
                  "TransmiR_Tissues_precursor", 
                `Transcription factors (TransmiR)` =
                  "TransmiR_Transcription_factors_precursor")
  
  ## test expected output
  expect_identical(available_categories, expected = expected)
  
})




