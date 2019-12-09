context("venn_fisher_test")

test_that("venn fisher test", {
  # monitor change
  mat.sig <- ezlimmaplot::ezvenn(mtt, p.cutoff=0.1, plot=FALSE)
  res <- venn_fisher_test(mat.sig)
  expect_equal(round(res["p-Value", "Overlapping of Up Regulation"], 5), 0.15283)
  expect_lte(res$'Overlapping of Up Regulation'[1], 0.153)
  expect_lte(res$'Overlapping of Up Regulation'[2], 9.727)
  expect_lte(res$'Overlapping of Up Regulation'[3], 0.302)
  expect_equal(res$'Overlapping of Up Regulation'[4], Inf)
  # no overlapping of down regulation in mat.sig
  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)

  mat.sig[,1] <- 1
  mat.sig[,2] <- 1
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Enrichment'[1], 1)
  expect_equal(res$'Overlapping of Enrichment'[2], 0)
  expect_equal(res$'Overlapping of Enrichment'[3], 0)
  expect_equal(res$'Overlapping of Enrichment'[4], Inf)

  mat.sig[,1] <- -1
  mat.sig[,2] <- 1
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Up Regulation'[1], 1)
  expect_equal(res$'Overlapping of Up Regulation'[2], 0)
  expect_equal(res$'Overlapping of Up Regulation'[3], 0)
  expect_equal(res$'Overlapping of Up Regulation'[4], Inf)
  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)

  mat.sig[,1] <- 1
  mat.sig[,2] <- 0
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Enrichment'[1], 1)
  expect_equal(res$'Overlapping of Enrichment'[2], 0)
  expect_equal(res$'Overlapping of Enrichment'[3], 0)
  expect_equal(res$'Overlapping of Enrichment'[4], Inf)

  mat.sig[,1] <- 0
  mat.sig[,2] <- -1
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Up Regulation'[1], 1)
  expect_equal(res$'Overlapping of Up Regulation'[2], 0)
  expect_equal(res$'Overlapping of Up Regulation'[3], 0)
  expect_equal(res$'Overlapping of Up Regulation'[4], Inf)
  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)

  mat.sig[,1] <- 0
  mat.sig[,2] <- 0
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Enrichment'[1], 1)
  expect_equal(res$'Overlapping of Enrichment'[2], 0)
  expect_equal(res$'Overlapping of Enrichment'[3], 0)
  expect_equal(res$'Overlapping of Enrichment'[4], Inf)

  mat.sig[,1] <- -1
  mat.sig[,2] <- -1
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Up Regulation'[1], 1)
  expect_equal(res$'Overlapping of Up Regulation'[2], 0)
  expect_equal(res$'Overlapping of Up Regulation'[3], 0)
  expect_equal(res$'Overlapping of Up Regulation'[4], Inf)
  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)
})
