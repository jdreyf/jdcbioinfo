context("z_2prop_pooled")

test_that("z_2prop_pooled", {

  x1 <- 3
  n1 <- 10
  x2 <- 8
  n2 <- 19

  expect_true(z_2prop_pooled(x1, x2, n1, n2, alternative="less") <
                z_2prop_pooled(x1, x2, n1, n2, alternative="greater"))

  expect_true(z_2prop_pooled(x1, x2, n1, n2, alternative="less", correct=FALSE) <
                z_2prop_pooled(x1, x2, n1, n2, alternative="greater", correct=FALSE))

  expect_equal(z_2prop_pooled(x1, x2, n1, n2, alternative="less")*2,
                z_2prop_pooled(x1, x2, n1, n2, alternative="two.sided"))

  expect_equal(z_2prop_pooled(x1, x2, n1, n2, alternative="less", correct=FALSE)*2,
               z_2prop_pooled(x1, x2, n1, n2, alternative="two.sided", correct=FALSE))

  expect_error(z_2prop_pooled(20, x2, n1, n2))

  expect_error(z_2prop_pooled(x1, x2, 1, n2))
})


