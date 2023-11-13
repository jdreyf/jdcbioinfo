test_that("add_avg_beta", {
  expect_silent(mtt.aab <- add_avg_beta(mtt))
  expect_error(add_avg_beta(mtt.aab))
  expect_gt(min(mtt.aab |> dplyr::select(tidyselect::ends_with("beta.avg"))), 0)
})
