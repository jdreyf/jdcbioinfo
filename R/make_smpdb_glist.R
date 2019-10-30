#' Make a feature list for SMPDB pathways.
#'
#' Make a feature list for SMPDB pathways.
#'
#'
#' @param smpdb.tab A dataframe of SMPDB data, whith column names in \code{id.cols.smpdb}, \code{pwys.name.col}, and \code{pwys.id.col}.
#' @param annot Metabolite annotation, with column names in \code{id.cols.annot}
#' @param id.cols.smpdb Column names in \code{smpdb.tab} for HMDB, ChEBI, and KEGG IDs.
#' @param id.cols.annot Column names in \code{annot} for HMDB, ChEBI, and KEGG IDs.
#' @param pwys.name.col Column name in \code{smpdb.tab} for pathway names.
#' @param pwys.id.col Column name in \code{smpdb.tab} for pathway IDs.
#' @param pwys.type.col Column name in \code{smpdb.tab} for pathway types.
#' @param pwys.type Type of pathways to use. Default is "Metabolic".
#' @export

make_smpdb_glist <- function(smpdb.tab, annot, id.cols.smpdb=c("Metabolite.ID", "ChEBI.ID", "KEGG.ID"), id.cols.annot=c("HMDB", "ChEBI", "KEGG"),
                             pwys.name.col="Pathway.Name", pwys.id.col="SMPDB.ID", pwys.type.col="Pathway.Type", pwys.type="Metabolic"){

  stopifnot(id.cols.smpdb %in% colnames(smpdb.tab), pwys.name.col %in% colnames(smpdb.tab), pwys.id.col %in% colnames(smpdb.tab),
            pwys.type.col %in% colnames(smpdb.tab))
  stopifnot(id.cols.annot %in% colnames(annot), !is.null(rownames(annot)))

  smpdb.tab <- smpdb.tab[, c(id.cols.smpdb, pwys.name.col, pwys.id.col, pwys.type.col)]
  smpdb.tab <- smpdb.tab[smpdb.tab[, pwys.type.col] %in% pwys.type]
  smpdb.tab[is.na(smpdb.tab) | smpdb.tab %in% c("NA", "N/A", "n/a")] <- ""

  if (mean(grepl("^HMBD",  smpdb.tab[,1])) < 0.5) {
    smpdb.tab[,1][smpdb.tab[,1]!=""] <- paste0("HMBD", smpdb.tab[,1][ smpdb.tab[,1]!=""])
  }
  if (mean(grepl("^CHEBI:",  smpdb.tab[,2])) < 0.5) {
    smpdb.tab[,2][smpdb.tab[,2]!=""] <- paste0("CHEBI:", smpdb.tab[,2][ smpdb.tab[,2]!=""])
  }
  if (mean(grepl("^C",  smpdb.tab[,3])) < 0.5) {
    smpdb.tab[,3][smpdb.tab[,3]!=""] <- paste0("C", smpdb.tab[,3][ smpdb.tab[,3]!=""])
  }

  annot <- annot[, id.cols.annot]
  pwys <- smpdb.tab[, pwys.name.col]
  G <- list()
  for (pwy in unique(pwys)) {
    smp <- unique(smpdb[smpdb.tab[, pwys.name.col] == pwy, pwys.id.col])
    stopifnot(length(smp)==1)

    hmdb.ids <- smpdb[pwys == pwy, 1]
    hmdb.ids <- hmdb.ids[hmdb.ids!=""]
    chebi.ids <- smpdb[pwys == pwy, 2]
    chebi.ids <- chebi.ids[chebi.ids != ""]
    kegg.ids <- smpdb[pwys == pwy, 3]
    kegg.ids <- kegg.ids[kegg.ids !=""]

    nm <- gsub(" |/", "_" , pwy)
    G[[nm]] <- list(name=nm, description=smp, genes=rownames(annot)[annot[,1] %in% hmdb.ids | annot[,2] %in% chebi.ids | annot[,3] %in% kegg.ids])
  }

 return(G)
}
