context("ezchip2chip")

test_that("ezchip2chip", {
  # want to convert "g1" to "gene1"
  G <- list(list(name="pwy1", description=NA, genes=paste0("g", 1:10)),
            list(name="pwy2", description=NA, genes=paste0("g", 11:20)),
            list(name="pwy3", description=NA, genes=paste0("g", 21:30)))

  annot.tmp <- data.frame(annot, sym=paste0("g", 1:100))
  rownames(annot.tmp) <- rownames(M)
  annot.tmp$sym[1] <- "g1 __ g11"

  ecc <- ezchip2chip(G=G, annot=annot.tmp, geneid.col = "sym", split.str = " __ ")

  expect_equal(length(ecc), 3)
  expect_equal(length(ecc[[1]][[3]]), 10)
  expect_equal(length(ecc[[2]][[3]]), 11)
  expect_equal(length(ecc[[3]][[3]]), 10)
  expect_true("gene1" %in% ecc[[1]][[3]])
  expect_true("gene1" %in% ecc[[2]][[3]])
  expect_false("gene1" %in% ecc[[3]][[3]])
})
