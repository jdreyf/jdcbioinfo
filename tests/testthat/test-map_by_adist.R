context("map_by_adist")

test_that("ignore case", {

  x <- c("ABC", "cde", "adcf")
  y <- c("cdf", "abc")

  tab1 <- map_by_adist(x=x, y=y, max.norm.dist=0.5)
  tab2 <- data.frame(x=c("cde", "ABC"), y=c("cdf", "abc"))
  expect_equal(tab1, tab2)
})

test_that("not ignore case", {

  x <- c("ABC", "cde", "adcf")
  y <- c("cdf", "abc")

  tab1 <- map_by_adist(x=x, y=y, ignore.case=FALSE, max.norm.dist=0.5)
  tab2 <- data.frame(x=c("cde", "adcf"), y=c("cdf", "abc"))
  expect_equal(tab1, tab2)
})
