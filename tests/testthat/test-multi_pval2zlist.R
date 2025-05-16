test_that("works with NAs", {
  mtt_na <- mtt
  mtt_na[1:2, c("Middle3vsFirst3.p", "Middle3vsFirst3.FDR", "Middle3vsFirst3.logFC", "Middle3vsFirst3.FC")] <- NA
  mtt_na[2, c("Last3vsFirst3.p", "Last3vsFirst3.FDR", "Last3vsFirst3.logFC", "Last3vsFirst3.FC")] <- NA
  z.lst <- multi_pval2zlist(mtt_na, sort = FALSE)
  expect_equal(sum(is.na(z.lst[[1]])), 0)
  expect_equal(sum(is.na(z.lst[[2]])), 0)
  expect_length(z.lst[[1]], 98)
  expect_length(z.lst[[2]], 99)

  z2 <- multi_pval2zlist(mtt_na[1:2,], sort = FALSE)
  expect_length(z2[[1]], 0)
  expect_length(z2[[2]], 1)
})
