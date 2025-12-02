test_that("multiple contrasts", {
  pheno$Group <- factor(pheno$Group)
  des <- model.matrix(~0 + Group, data = pheno)
  colnames(des) <- gsub("Group", "", colnames(des))
  dds <- DESeq2::DESeqDataSetFromMatrix(round(dge$counts), colData = pheno, design = des)

  res <- deseq2_contrasts(dds, grp=pheno$Group, contrast.v = contr.v, shrunken = FALSE, parallel = FALSE)
  expect_true(all(rownames(res)[1:3] %in% paste0("gene", 1:3)))

  res2 <- deseq2_contrasts(dds, grp=pheno$Group, contrast.v = contr.v, shrunken = TRUE, parallel=FALSE)
  expect_true(all(rownames(res2)[1:3] %in% paste0("gene", 1:3)))

  res3 <- deseq2_contrasts(dds, grp=pheno$Group, contrast.v = contr.v, independentFiltering = TRUE, add_cooks=TRUE, parallel = FALSE)
  expect_true(all(rownames(res3)[1:3] %in% paste0("gene", 1:3)))
})
