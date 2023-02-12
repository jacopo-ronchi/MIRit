## code to prepare `geneCounts` and `mirnaCounts` datasets goes here

## downlaod and extract supplementary files from GEO
GEOquery::getGEOSuppFiles(GEO = "GSE63511")
untar("GSE63511/GSE63511_RAW.tar", exdir = "GSE63511")

## list files in data directory
rf <- list.files("GSE63511", full.names = TRUE)

## list RNA files and miRNA files
rf_g <- rf[grep(".RNA", rf, fixed = TRUE)]
rf_m <- rf[grep("smRNA", rf, fixed = TRUE)]

## retain only miRNA files with matched RNA files
rf_m <- rf_m[1:16]

## read counts to form count matrices for RNA and miRNAs
g <- sapply(rf_g, function(x) {
  a <- read.csv(x, sep = "\t", header = FALSE)
  a$V2
})
rownames(g) <- read.csv(rf_g[1], sep = "\t", header = FALSE)$V1

m <- sapply(rf_m, function(x) {
  a <- read.csv(x, sep = "\t", header = FALSE)
  a$V2
})
rownames(m) <- read.csv(rf_m[1], sep = "\t", header = FALSE)$V1

## remove ambiguous last 5 lines
g <- g[-c((nrow(g) - 4):nrow(g)), ]
m <- m[-c((nrow(m) - 4):nrow(m)), ]

## name samples and create metadata data.frame
meta <- data.frame(id = character(16), disease = character(16))
meta$id <- c(paste("PTC", seq(8)), paste("NTH", seq(8)))
meta$disease <- factor(c(rep("PTC", 8), rep("NTH", 8)))
meta$patient <- factor(c(rep(paste("patient_", seq(8), sep = ""), 2)))
colnames(m) <- colnames(g) <- meta$id

## remove the downloaded GEO data folder
unlink("GSE63511", recursive = TRUE)

## rename datasets
geneCounts <- g
mirnaCounts <- m

## create package data
usethis::use_data(geneCounts, overwrite = TRUE)
usethis::use_data(mirnaCounts, overwrite = TRUE)
