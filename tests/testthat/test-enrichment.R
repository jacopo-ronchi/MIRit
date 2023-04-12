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



test_that("MicroRNA ORA works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform ORA of DE-miRNAs
  enr <- enrichMirnas(obj,
                      organism = "Homo sapiens",
                      category = "GO_Annotations_mature",
                      pCutoff = 0.1,
                      pAdjustment = "fdr",
                      minHits = 3,
                      mirnaType = "mature")
  
  ## test enrichment results
  expect_snapshot(enrichmentResults(enr$upregulated))
  
})



test_that("MicroRNA GSEA works", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform GSEA of miRNAs
  gse <- gseaMirnas(obj,
                    organism = "Homo sapiens",
                    category = "GO_Annotations_mature",
                    pCutoff = 0.1,
                    pAdjustment = "fdr",
                    minHits = 3,
                    mirnaType = "mature")
  
  ## test enrichment results
  expect_snapshot(enrichmentResults(gse))
  
})



test_that("Intergated targets ORA works (with GO)", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform ORA of integrated targets
  enr <- enrichTargets(obj,
                       integratedTargets = TRUE,
                       database = "GO",
                       organism = "Homo sapiens",
                       simplifyGO = FALSE)
  
  ## test enrichment results
  expect_snapshot(enrichmentResults(enr$`UP-miRNA targets`))
  
})



test_that("Genes ORA works (with Reactome, DisGeNet and WikiPathways)", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform ORA - Reactome of genes
  rc <- enrichGenes(obj,
                    database = "Reactome",
                    organism = "Homo sapiens")
  
  ## test ORA - Reactome results
  expect_snapshot(enrichmentResults(rc$downregulated))
  
  ## perform ORA - DisGeNet of genes
  dg <- enrichGenes(obj,
                    database = "KEGG",
                    organism = "Homo sapiens",
                    pCutoff = 0.05,
                    pAdjustment = "none")
  
  ## test ORA - DisGeNet results
  expect_snapshot(enrichmentResults(dg$downregulated))
  
  ## perform ORA - WikiPathways of genes
  wp <- enrichGenes(obj,
                    database = "WikiPathways",
                    organism = "Homo sapiens")
  
  ## test ORA - WikiPathways results
  expect_snapshot(enrichmentResults(wp$downregulated))
  
})



test_that("Genes GSEA works (with KEGG)", {
  
  ## load the example MirnaExperiment object
  obj <- loadExamples()
  
  ## perform GSEA - KEGG of genes
  kg <- gseaGenes(obj,
                  database = "KEGG",
                  organism = "Homo sapiens")
  
  ## test GSEA - KEGG results
  expect_snapshot(enrichmentResults(kg))
  
})
