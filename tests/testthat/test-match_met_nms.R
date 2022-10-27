context("match metabolite names")

test_that("match_met_nms", {
  x <- c("a", "b", "cate", "2 Methyl-malonyl CoA", "adate", "adic acid")
  tble <- c("aa", "a",  "cic acid", "adate", "2 methylmalonyl coa")
  expect_equal(match_met_nms(x, tble), c(2, NA, 3, 5, 4, 4))
})


