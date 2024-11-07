context("edger_contrasts")

test_that("edger contrasts", {
  #default
  res <- edger_contrasts(dge, grp=grp, contrast.v = contr.v)
  expect_equal(rownames(res$mtt)[1:3], paste0("gene", c(3,2,1)))

  #LRT
  suppressWarnings(res <- edger_contrasts(dge, grp=grp, contrast.v = contr.v, test="LRT"))
  expect_equal(rownames(res$mtt)[1:3], paste0("gene", c(3,62,1)))

  # design
  res <- edger_contrasts(dge, grp=grp, contrast.v = contr.v, design = design)
  expect_equal(rownames(res$mtt)[1:3], paste0("gene", 3:1))

  # weights
  dge1 <- dge
  dge1$weights <- el$weights
  suppressWarnings(res <- edger_contrasts(dge1, grp=grp, contrast.v = contr.v, test="WeightedFT"))
  expect_equal(rownames(res$mtt)[1:3], paste0("gene", 1:3))
})
