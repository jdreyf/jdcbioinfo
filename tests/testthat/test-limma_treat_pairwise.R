context("limma_treat_pairwise")

test_that("Matrix as input", {
  #default
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp")
  expect_true(all(rownames(ttf)[1:3] %in% paste0("gene", 1:3)))

  #trend
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp", trend = TRUE)
  expect_true(all(rownames(ttf)[1:3] %in% paste0("gene", 1:3)))

  # design
  ttf <- limma_treat_pairwise(M, grp=grp, prefix = "grp", design = design)
  expect_true(all(rownames(ttf)[1:3] %in% paste0("gene", 1:3)))
})

test_that("DGEList as input", {
  #default
  ttf <- limma_treat_pairwise(el, grp=grp, prefix = "grp")
  expect_true(all(rownames(ttf)[1:3] %in% paste0("gene", 1:3)))

  # design
  ttf <- limma_treat_pairwise(el, grp=grp, prefix = "grp", design = design)
  expect_true(all(rownames(ttf)[1:3] %in% paste0("gene", 1:3)))
})
