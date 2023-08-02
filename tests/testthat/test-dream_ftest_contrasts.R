context("dream_ftest_contrasts")

test_that("test dream_ftest_contrasts", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ 0 + Batch + (1|Individual) + (1|Tissue)
  L1 <- getContrast( geneExpr, form, info, c("BatchB1", "BatchB2"))
  L2 <- getContrast( geneExpr, form, info, c("BatchB1", "BatchB3"))
  L3 <- getContrast( geneExpr, form, info, c("BatchB2", "BatchB3"))
  L <- cbind(L1, L2, L3)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  fit1 <- eBayes(fit1)
  tt1 <- topTable(fit1, coef=c("L1", "L2", "L3"), sort.by="none")
  tt1 <- tt1[order(tt1$P.Value), ]

  contr.v <- c(Batch1vs2 = "BatchB1 - BatchB2", Batch1vs3 = "BatchB1 - BatchB3", Batch2vs3 = "BatchB2 - BatchB3")
  tt2 <- dream_ftest_contrasts(geneExpr[1:10,], formula=form, pheno=info, contrast.v=contr.v, grp=info$Batch, prefix = "Group")
  expect_equal(signif(tt1$P.Value, 4), signif(tt2$Group.p, 4))
})

