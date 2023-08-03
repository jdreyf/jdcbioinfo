context("dream_contrasts")

test_that("compare to dream", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ Batch + (1|Individual) + (1|Tissue)
  contr.v <- c(Batch2vs3 = "BatchB2 - BatchB3")
  L <- makeContrastsDream(form, data=info, contrasts=contr.v)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  fit1 <- eBayes(fit1)
  tt1 <- topTable(fit1, coef="Batch2vs3", sort.by="p")

  tt2 <- dream_contrasts(geneExpr[1:10,], formula=form, pheno=info, contrast.v=contr.v, grp=info$Batch)
  expect_equal(tt1$P.Value, tt2$Batch2vs3.p)
})

