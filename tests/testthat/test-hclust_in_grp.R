context("hclust_in_grp")

test_that("centered", {
  grps <- rep(LETTERS[1:4], each = 25)
  lbs <- hclust_in_grp(M, grp = grps, sc="ctr")
  expect_equal(as.vector(lbs[1:3]), c("gene3", "gene1", "gene2"))
})

test_that("z-scored", {
  grps <- rep(LETTERS[1:4], each = 25)
  lbs <- hclust_in_grp(M, grp = grps, sc="z")
  expect_equal(as.vector(lbs[1:3]), c("gene16", "gene21", "gene2"))
})
