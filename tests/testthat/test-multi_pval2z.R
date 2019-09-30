context("multi_pval2z")

mtt <- ezlimma::limma_contrasts(M, grp=grp, contrast.v = contr.v)
mtt1 <- mtt
mtt1$Middle3vsFirst3.Direction <- ifelse(mtt1$Middle3vsFirst3.logFC>0, "Up", "Down")
mtt1$Last3vsFirst3.Direction <- ifelse(mtt1$Last3vsFirst3.logFC>0, "Up", "Down")

test_that("only tab", {
  res <- multi_pval2z(mtt)
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})

test_that("with prefix.v", {
  res <- multi_pval2z(mtt, prefix.v = names(contr.v))
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})

test_that("one comparison", {
  res <- multi_pval2z(mtt, prefix.v = names(contr.v)[1])
  expect_equal(as.vector(res[1:3]), c(7.972129881, 1.426847696, 6.473577613))
})

test_that("with Up/Down", {
  res <- multi_pval2z(mtt1, direction.suffix = "Direction")
  expect_equal(as.vector(res[1, ]), c(7.972129881, 7.96543488))
})
