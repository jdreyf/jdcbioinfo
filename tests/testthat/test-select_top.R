context("select_top")

test_that("without prefix.v", {
  topg <- select_top(mtt, ntop=3)
  expect_equal(topg, c("gene3", "gene1", "gene68", "gene2", "gene5" ))
})

test_that("with prefix.v", {
  topg <- select_top(mtt, prefix.v=names(contr.v), ntop=3)
  expect_equal(topg, c("gene3", "gene1", "gene68", "gene2", "gene5" ))
})

test_that("not each", {
  topg <- select_top(mtt, ntop=3, each = FALSE)
  expect_equal(topg, c("gene3", "gene2", "gene1"))
})

test_that("large ntop", {
  topg <- select_top(mtt, ntop=200)
  expect_equal(topg, rownames(mtt))
})
