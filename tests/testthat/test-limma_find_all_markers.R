test_that("multiplication works", {
  lfam <- limma_find_all_markers(object = M, grp = grp)
  expect_equal(nrow(lfam), nrow(M))
  expect_equal(ncol(lfam), 3 + 6*length(unique(grp)))
  expect_lte(lfam["gene1", "Middle3.up.p"], 0.01)
  expect_lte(lfam["gene2", "Last3.up.p"], 0.01)
  expect_lte(lfam["gene3", "First3.down.p"], 0.01)
})
