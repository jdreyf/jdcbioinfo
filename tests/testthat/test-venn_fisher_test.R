context("venn_fisher_test")

test_that("venn fisher test", {
  mat.sig <- ezlimmaplot::ezvenn(mtt, p.cutoff=0.1, plot=FALSE)
  sig.up <- mat.sig
  sig.up[sig.up != 1] <- 2
  sig.up <- as.data.frame(sig.up)
  sig.up[,1] <- factor(sig.up[,1], levels=c(1,2))
  sig.up[,2] <- factor(sig.up[,2], levels=c(1,2))
  tab.up <- table(sig.up[,1], sig.up[,2])
  res.up <- stats::fisher.test(tab.up, alternative="greater")
  res <- venn_fisher_test(mat.sig, directional = TRUE)
  # test for matrix of (-1, 0, 1)
  expect_equal(res$'Overlapping of Up Regulation'[1], res.up$p.value)
  expect_equal(res$'Overlapping of Up Regulation'[2], as.numeric(res.up$estimate))
  expect_lte(res$'Overlapping of Up Regulation'[3], res.up$conf.int[1])
  expect_equal(res$'Overlapping of Up Regulation'[4], res.up$conf.int[2])

  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)

  # test for matrix of (-1, 0, 1) with  directional = FALSE
  expect_error(venn_fisher_test(mat.sig, directional = FALSE))

  # test for matrix of (0, 1) with  directional = TRUE/FALSE
  mat.sig[,1] <- 1
  mat.sig[,2] <- 0
  res <- venn_fisher_test(mat.sig)
  expect_equal(res$'Overlapping of Up Regulation'[1], 1)
  expect_equal(res$'Overlapping of Up Regulation'[2], 0)
  expect_equal(res$'Overlapping of Up Regulation'[3], 0)
  expect_equal(res$'Overlapping of Up Regulation'[4], Inf)
  expect_equal(res$'Overlapping of Down Regulation'[1], 1)
  expect_equal(res$'Overlapping of Down Regulation'[2], 0)
  expect_equal(res$'Overlapping of Down Regulation'[3], 0)
  expect_equal(res$'Overlapping of Down Regulation'[4], Inf)

  res <- venn_fisher_test(mat.sig, directional = FALSE)
  expect_equal(res$'Overlapping of Enrichment'[1], 1)
  expect_equal(res$'Overlapping of Enrichment'[2], 0)
  expect_equal(res$'Overlapping of Enrichment'[3], 0)
  expect_equal(res$'Overlapping of Enrichment'[4], Inf)

  # test for matrix of (0, 0) with  directional = TRUE/FALSE
  mat.sig[,1] <- 0
  mat.sig[,2] <- 0
  expect_error(venn_fisher_test(mat.sig, directional = FALSE))
  expect_error(venn_fisher_test(mat.sig, directional = TRUE))
})
