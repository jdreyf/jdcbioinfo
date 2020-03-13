context("limma_treat_max")

test_that("Matrix as input", {

  #default
  ttf <- limma_treat_max(M, grp=grp, contrast.v = contr.v, prefix = "grp")
  expect_equal(rownames(ttf)[1:3], paste0("gene", 3:1))

  #trend
  ttf <- limma_treat_max(M, grp=grp, contrast.v = contr.v, prefix = "grp", trend = TRUE)
  expect_equal(rownames(ttf)[1:3], paste0("gene", c(3,1,2)))

  # design
  ttf <- limma_treat_max(M, grp=grp, contrast.v = contr.v, prefix = "grp", design = design)
  expect_equal(rownames(ttf)[1:3], paste0("gene", 3:1))

  # t-test
  ttf <- limma_treat_max(M, grp=grp, contrast.v = contr.v[1], prefix = "grp")
  expect_equal(rownames(ttf)[1:2], paste0("gene", c(3,1)))
})

test_that("DGEList as input", {

  #default
  ttf <- limma_treat_max(el, grp=grp, contrast.v = contr.v, prefix = "grp")
  expect_equal(rownames(ttf)[1:3], paste0("gene", c(3,1,2)))

  # design
  ttf <- limma_treat_max(el, grp=grp, contrast.v = contr.v, prefix = "grp", design = design)
  expect_equal(rownames(ttf)[1:3], paste0("gene",  c(3,1,2)))

  # t-test
  ttf <- limma_treat_max(el, grp=grp, contrast.v = contr.v[1], prefix = "grp")
  expect_equal(rownames(ttf)[1:2], paste0("gene", c(3,1)))
})
