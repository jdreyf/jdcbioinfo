context("hclust_and_heat")

test_that("default", {
  clus_df <- hclust_and_heat(M, annot=annot, pheno.df = pheno, plot=FALSE, verbose=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("gene79", "gene13", "gene53"))
  expect_equal(length(unique(clus_df$Cluster)), 6)
})


test_that("deepSplit=2", {
  clus_df <- hclust_and_heat(M, annot=annot, pheno.df = pheno, deepSplit=2,minClusterSize=5,  plot=FALSE, verbose=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("gene79", "gene13", "gene53"))
  expect_equal(length(unique(clus_df$Cluster)), 12)
})


test_that("clip=1", {
  clus_df <- hclust_and_heat(M, annot=annot, pheno.df = pheno, clip=1, plot=FALSE, verbose=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("gene49", "gene71", "gene25"))
  expect_equal(length(unique(clus_df$Cluster)), 5)
})
