context("binary_str_lengths")

test_that("Examples", {
  expect_equal(length(binary_str_lengths(0)), 0)
  expect_equal(binary_str_lengths(1), 1)
  expect_equal(binary_str_lengths(c(1, 1)), 2)
  expect_equal(binary_str_lengths(c(1, 1, 1, 0, 0, 1, 0)), c(1, 3))
  expect_error(binary_str_lengths(TRUE))
  expect_error(binary_str_lengths(c("1")))
})
