#' Get gene anotation from biomart
#'
#' Get gene anotation from biomart.
#'
#' @param ids Ensembl gene IDs
#' @param organism "hsa" or "mmu" for human or mouse
#' @return Data frame with gene anotation.
#' @export
#'
ezgetBM <- function(ids, organism=c("hsa", "mmu")){

  organism <- match.arg(organism)

  if (organism=="mmu") {
    ensembl <- biomaRt::useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
    annot <- biomaRt::getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "entrezgene_id", "description", "gene_biotype"),
                  filters="ensembl_gene_id", values=ids, mart=ensembl)
  } else if (organism=="hsa") {
    ensembl <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    annot <- biomaRt::getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "description", "gene_biotype"),
                  filters="ensembl_gene_id", values=ids, mart=ensembl)
  }

  #need to collapse multi ensembl ids
  colnames(annot) <- c("ensembl", "symbol", "entrez", "description", "gene_biotype")
  annot$description <- gsub(" \\[Source.+", "",  annot$description)
  annot <- stats::aggregate(.~ensembl, data = annot, FUN = function(x) paste(sort(unique(x)), collapse = " /// "))
  annot <- data.frame(annot, row.names = "ensembl")

  return(annot)
}
