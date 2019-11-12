context("signed_rankprod")

test_that("same direction", {
  mat_z <- multi_pval2z(mtt)
  rp <- signed_rankprod(mat_z, nsim = 1e5-1, same.dirction = TRUE)
  expect_equal(rownames(rp)[1:3], c("gene3", "gene41", "gene2"))
})

test_that("opposite direction", {
  mat_z <- multi_pval2z(mtt)
  rp <- signed_rankprod(mat_z, nsim = 1e5-1, same.dirction = FALSE)
  expect_equal(rownames(rp)[1:3], c("gene3", "gene41", "gene2"))
})
