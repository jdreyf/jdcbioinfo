library(covr)
library(ppcor)
library(limma)
library(ezlimma)
library(ezlimmaplot)
library(testthat)
library(variancePartition)
library(BiocParallel)

data(varPartData)

suppressWarnings(set.seed(100, sample.kind = "Rounding"))
means <- sample(seq(1, 20, by = 0.01), 100, replace = TRUE)
M <-t(sapply(1:100, FUN = function(i) rnorm(n = 9, mean = means[i], sd = 1/means[i])))
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
M[1, 4:6] <- M[1, 4:6] + mean(M[1,])
M[2, 7:9] <- M[2, 7:9] + mean(M[2,])
M[3, 4:9] <- M[3, 4:9] + mean(M[3,])

covar1 <- sample(0:1, 9, replace = TRUE)
covar2 <- rnorm(ncol(M))
design <- cbind(First3=c(1,1,1,0,0,0,0,0,0), Middle3 = c(0,0,0,1,1,1,0,0,0),
                Last3=c(0,0,0,0,0,0,1,1,1), Covar1 = covar1, Covar2 = covar2)

dge <- edgeR::DGEList(counts = 2^M)
dge <- edgeR::calcNormFactors(dge)
el <- limma::voom(dge, design=design, plot = FALSE)

grp <- rep(c("First3", "Middle3", "Last3"), each=3)
contr.v <- c(Middle3vsFirst3="Middle3-First3", Last3vsFirst3="Last3-First3")

pval <- runif(100)
logfc <- sample(seq(-2, 2, by = 0.01), 100, replace = TRUE)
direction <- sample(c("Up", "Down"), 100, replace = TRUE)

mtt <- ezlimma::limma_contrasts(M, grp=grp, contrast.v = contr.v)

annot <- data.frame(Gene=rownames(M), row.names=rownames(M))
pheno <- data.frame(Group=grp, row.names=colnames(M))
