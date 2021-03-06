% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ezgetBM.R
\name{ezgetBM}
\alias{ezgetBM}
\title{Get gene anotation from biomart}
\usage{
ezgetBM(ids, organism = c("hsa", "mmu"))
}
\arguments{
\item{ids}{Ensembl gene IDs}

\item{organism}{"hsa" or "mmu" for human or mouse}
}
\value{
Data frame with gene anotation.
}
\description{
Get gene anotation from biomart.
}
