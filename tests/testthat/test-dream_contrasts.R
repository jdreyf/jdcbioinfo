context("dream_contrasts")

test_that("test dream_contrasts", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ Batch + (1|Individual) + (1|Tissue)
  L <- getContrast( geneExpr, form, info, c("BatchB2", "BatchB3"))
  set.seed(100)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  tt1 <- topTable(fit1, coef = "L1", sort.by = "p")

  contr.v <- c(Batch2vs3 = "BatchB2 - BatchB3")
  set.seed(100)
  tt2 <- dream_contrasts(geneExpr[1:10,], formula=form, pheno=info, contrast.v = contr.v, grp = info$Batch)
  expect_equal(tt1$P.Value, tt2$Batch2vs3.p)
})

