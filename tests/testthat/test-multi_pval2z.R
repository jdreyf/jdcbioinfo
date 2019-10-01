context("multi_pval2z")

test_that("without prefix.v", {
  res <- multi_pval2z(mtt)
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})

test_that("with prefix.v", {
  res <- multi_pval2z(mtt, prefix.v = names(contr.v))
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})

test_that("Invalid p-value input", {
  mtt_invalid <- mtt
  mtt_invalid$Last3vsFirst3.p[1] <- 2
  expect_error(multi_pval2z(mtt_invalid))
})

test_that("Direction input", {
  mtt_d <- mtt
  mtt_d$Middle3vsFirst3.Direction <- ifelse(mtt_d$Middle3vsFirst3.logFC>0, "Up", "Down")
  mtt_d$Last3vsFirst3.Direction <- ifelse(mtt_d$Last3vsFirst3.logFC>0, "Up", "Down")
  res <- multi_pval2z(mtt_d, direction.suffix="Direction")
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})
