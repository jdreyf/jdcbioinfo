context("ezbacon")

test_that("small example", {
  mtt.bc <- ezlimma::limma_contrasts(M, grp=grp, contrast.v = contr.v, cols=c("logFC", "z"))
  mtt.bc <- jdcbioinfo::ezbacon(tab=mtt.bc)
  expect_gte(cor(mtt.bc$Middle3vsFirst3.z, mtt.bc$Middle3vsFirst3_bacon.z), 0.99)
  expect_gte(cor(mtt.bc$Last3vsFirst3.z, mtt.bc$Last3vsFirst3_bacon.z), 0.99)
})
