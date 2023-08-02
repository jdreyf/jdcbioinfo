context("dream_ftest_pairwise")

test_that("test dream_ftest_pairwise", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ 0 + Batch + (1|Individual) + (1|Tissue)
  L1 <- getContrast( geneExpr, form, info, c("BatchB1", "BatchB2"))
  L2 <- getContrast( geneExpr, form, info, c("BatchB1", "BatchB3"))
  L3 <- getContrast( geneExpr, form, info, c("BatchB1", "BatchB4"))
  L4 <- getContrast( geneExpr, form, info, c("BatchB2", "BatchB3"))
  L5 <- getContrast( geneExpr, form, info, c("BatchB2", "BatchB4"))
  L6 <- getContrast( geneExpr, form, info, c("BatchB3", "BatchB4"))
  L <- cbind(L1, L2, L3, L4, L5, L6)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  fit1 <- eBayes(fit1)
  tt1 <- topTable(fit1, coef=colnames(L), sort.by="none")
  tt1 <- tt1[order(tt1$P.Value), ]

  tt2 <- dream_ftest_pairwise(geneExpr[1:10,], formula=form, pheno=info, grp=paste0("Batch", info$Batch), prefix = "Group")
  expect_equal(signif(tt1$P.Value, 4), signif(tt2$Group.p, 4))
})

