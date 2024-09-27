context("limma_find_markers")

test_that("default", {
  res <- suppressWarnings(limma_find_markers(M, grp=grp, nsim=100))
  expect_equal(rownames(res)[which.min(res$Last3.up.p)], "gene2")
  expect_equal(rownames(res)[which.min(res$Middle3.up.p)], "gene1")
})

test_that("down", {
  res <- suppressWarnings(limma_find_markers(M, grp=grp, direction = "down", nsim=100))
  expect_equal(rownames(res)[which.min(res$First3.down.p)], "gene3")
  expect_equal(rownames(res)[which.min(res$Middle3.down.p)], "gene18")
})
