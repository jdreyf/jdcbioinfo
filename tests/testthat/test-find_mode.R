context("find_mode")

test_that("find mode", {
  set.seed(100)
  x1 <- stats::rnorm(1000, mean=10, sd=2)
  x2 <- stats::rnorm(1000, mean=20, sd=3)
  x <- c(x1, x2)
  mds <- find_mode(x, bin=30, plot=FALSE)
  expect_equal(round(mds$mds.vals), c(10, 20))
})


