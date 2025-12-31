test_that("add_avg_beta", {
  expect_silent(mtt.aab <- add_avg_beta(mtt))
  expect_equal(mtt.aab[1,1], 2^mtt[1,1]/(1 + 2^mtt[1,1])) # don't want to import sesame
  expect_error(add_avg_beta(mtt.aab))
  expect_gt(min(mtt.aab |> dplyr::select(tidyselect::ends_with("beta.avg"))), 0)
})
