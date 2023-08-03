context("dream_ftest_pairwise")

test_that("test dream_ftest_pairwise", {
  info$Batch <- paste0("B", info$Batch)
  form <- ~ 0 + Batch + (1|Individual) + (1|Tissue)
  grp <- paste0("Batch", info$Batch)
  comb <- utils::combn(unique(grp), 2)
  contr.v <- character(0)
  for(i in 1:ncol(comb)){
    contr.v[paste0(comb[2, i], "vs", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }
  L <- makeContrastsDream(form, data=info, contrasts=contr.v)
  fit1 <- dream(geneExpr[1:10,], formula=form, data=info, L=L)
  fit1 <- eBayes(fit1)
  tt1 <- topTable(fit1, coef=colnames(L), sort.by="none")
  tt1 <- tt1[order(tt1$P.Value), ]

  tt2 <- dream_ftest_pairwise(geneExpr[1:10,], formula=form, pheno=info, grp=grp, prefix = "Group")
  expect_equal(tt1$P.Value, tt2$Group.p)
})

