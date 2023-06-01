context("hclust_and_heat")

test_that("clip", {
  clus_df <- hclust_and_heat(M, annot=annot, pheno.df = pheno, plot=FALSE, verbose=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("gene79", "gene13", "gene53"))
  expect_equal(length(unique(clus_df$Cluster)), 6)

  # clip after clustering
  clus_df2 <- hclust_and_heat(M, annot=annot, pheno.df = pheno, plot=FALSE, verbose=FALSE, clip=1)
  expect_equal(clus_df2[rownames(clus_df), "Cluster"], clus_df$Cluster)
  expect_equal(length(unique(clus_df$Cluster)), length(unique(clus_df2$Cluster)))
})


test_that("deepSplit=2", {
  clus_df <- hclust_and_heat(M, annot=annot, pheno.df = pheno, deepSplit=2,minClusterSize=5,  plot=FALSE, verbose=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("gene79", "gene13", "gene53"))
  expect_equal(length(unique(clus_df$Cluster)), 12)
})

