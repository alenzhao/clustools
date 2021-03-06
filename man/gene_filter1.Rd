% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/clustools.R
\name{gene_filter1}
\alias{gene_filter1}
\title{First filtering step}
\usage{
gene_filter1(d, min.cells, max.cells, min.reads)
}
\arguments{
\item{d}{Expression matrix with rows as genes and columns as cells}

\item{min.cells}{Minimum number of cells in which a given gene is expressed}

\item{max.cells}{Maximum number of cells in which a given gene is expressed}

\item{min.reads}{Minimum number of reads per gene per cell}
}
\value{
Filtered expression matrix in which only genes that are expressed in
more than \code{min.cells} with more than \code{min.reads} reads and also are
expressed in less than [total number of cells - \code{max.cells}].
}
\description{
Filter genes that would not contribute to clustering, because they are either
expreseed or not expressed in almost all cells.
}
\examples{
gene_filter1(quake, 3, 3, 2)
}

