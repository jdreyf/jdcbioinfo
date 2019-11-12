context("venn_fisher_test")

test_that("venn fisher test", {
  mat.sig <- ezlimmaplot::ezvenn(mtt, p.cutoff=0.1, plot=FALSE)
  res <- venn_fisher_test(mat.sig)
  expect_equal(round(res["p-Value", "Overlapping of Up Regulation"], 5), 0.15283)
})
