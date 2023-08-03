context("dream_ftest_contrasts")

test_that("test dream_ftest_contrasts", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ 0 + Batch + (1|Individual) + (1|Tissue)
  contr.v <- c(Batch1vs2 = "BatchB1 - BatchB2", Batch1vs3 = "BatchB1 - BatchB3", Batch2vs3 = "BatchB2 - BatchB3")
  L <- makeContrastsDream(form, data=info, contrasts=contr.v)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  fit1 <- eBayes(fit1)
  tt1 <- topTable(fit1, coef=colnames(L), sort.by="none")
  tt1 <- tt1[order(tt1$P.Value), ]

  tt2 <- dream_ftest_contrasts(geneExpr[1:10,], formula=form, pheno=info, contrast.v=contr.v, grp=info$Batch, prefix = "Group")
  expect_equal(tt1$P.Value, tt2$Group.p)
})

