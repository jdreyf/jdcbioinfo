context("limma_treat_pairwise")

test_that("Matrix as input", {
  #default
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp")
  expect_equal(rownames(ttf)[1:3], paste0("gene", 3:1))

  #trend
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp", trend = TRUE)
  expect_equal(rownames(ttf)[1:3], paste0("gene", c(3,1,2)))

  # design
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp", design = design)
  expect_equal(rownames(ttf)[1:3], paste0("gene", 3:1))
})

test_that("DGEList as input", {
  #default
  ttf <- limma_treat_pairwise(el, grp=grp, prefix = "grp")
  expect_equal(rownames(ttf)[1:3], paste0("gene", c(3,1,2)))

  # design
  ttf <- limma_treat_pairwise(el, grp=grp, prefix = "grp", design = design)
  expect_equal(rownames(ttf)[1:3], paste0("gene", c(3,1,2)))
})
