context("correct_null_distribution")

test_that("without prefix.v", {
  suppressWarnings(mtt1 <- correct_null_distribution(mtt))
  expect_equal(mtt1[1, "Middle3vsFirst3.p"], 1.559628e-15)
})

test_that("with prefix.v", {
  suppressWarnings(mtt1 <- correct_null_distribution(mtt, prefix.v = names(contr.v)))
  expect_equal(mtt1[1, "Middle3vsFirst3.p"], 1.559628e-15)
})


test_that("Direction input", {
  mtt_d <- mtt
  mtt_d$Middle3vsFirst3.Direction <- ifelse(mtt_d$Middle3vsFirst3.logFC>0, "Up", "Down")
  mtt_d$Last3vsFirst3.Direction <- ifelse(mtt_d$Last3vsFirst3.logFC>0, "Up", "Down")
  suppressWarnings(mtt_d <- correct_null_distribution(mtt_d, direction.suffix="Direction"))
  expect_equal(mtt_d[1, "Middle3vsFirst3.p"], 1.559628e-15)
})