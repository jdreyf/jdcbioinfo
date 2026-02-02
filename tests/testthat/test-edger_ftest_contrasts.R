context("edger_ftest_contrasts")

test_that("edger ftest contrasts", {

  #default
  res <- edger_ftest_contrasts(dge, grp=grp, contrast.v = contr.v)
  expect_equal(rownames(res$mtt)[1:4], paste0("gene", c(3,2,1,18)))

  #LRT
  res <- edger_ftest_contrasts(dge, grp=grp, contrast.v = contr.v, test="LRT")
  expect_equal(rownames(res$mtt)[1:4], paste0("gene", c(3,1,2,18)))

  # design
  res <- edger_ftest_contrasts(dge, grp=grp, contrast.v = contr.v, design = design)
  expect_equal(rownames(res$mtt)[1:4], paste0("gene", c(3,2,1,22)))

  # weights
  dge1 <- dge
  dge1$weights <- el$weights
  res <- edger_ftest_contrasts(dge1, grp=grp, contrast.v = contr.v, test="WeightedFT")
  # expect_equal(rownames(res$mtt)[1:4], paste0("gene", c(3,1,2,18)))
  # changed in edgR 4.0
  expect_equal(rownames(res$mtt)[1:4], paste0("gene", c(3,1,2,85)))
})
