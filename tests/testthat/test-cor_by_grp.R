context("cor_by_grp")

test_that("default", {
  ct <- cor_by_grp(M, phenotype=covar2, grp=grp)
  expect_equal(ct$Combined.p[1:3], c(0.02905244, 0.02926561, 0.03701487))
  expect_equal(rownames(ct)[1], "gene74")
})

test_that("spearman", {
  ct <- cor_by_grp(M, phenotype=covar2, grp=grp, method="spearman")
  expect_equal(rownames(ct)[1], "gene7")
})
