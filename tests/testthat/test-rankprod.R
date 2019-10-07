context("rankprod")

test_that("small example", {
  mm <- matrix(1:9, nrow=3)
  expect_error(rankprod(mm))
  dimnames(mm) <- list(letters[1:nrow(mm)], LETTERS[1:nrow(mm)])
  set.seed(1)
  rp <- rankprod(mm, nsim=6)
  expect_equal(rownames(rp), c("c", "b", "a"))
})
