% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reference.R
\name{generateReference}
\alias{generateReference}
\title{Generate reference data}
\usage{
generateReference(
  singleCellExpr,
  singleCellLabels,
  types = c("markers", "sigGenes", "signature", "cellTypeExpr"),
  log2Threshold = 1
)
}
\arguments{
\item{singleCellExpr}{Single-cell expression matrix in non-log scale. Rows are genes and columns are cells.}

\item{singleCellLabels}{An array contains cell type label for each cell. The order of labels must be the same as the order of cells in in singleCellExpr.}

\item{types}{An array of different type of reference to generate.}

\item{log2Threshold}{The threshold for log2FC for filtering markers.}
}
\value{
A list with names/keys are the same as in types.
If type is `markers`, the value is a list with names/keys are unique cell types and values for each key are an array of markers for the corresponding cell type.
If type is `sigGenes`, the value is the unique concatenation of all `markers`.
If type is `signature`, the value is an matrix of `sigGenes` x `cellTypes`.
If type is `cellTypeExpr`,  the value is an matrix of `all genes` x `cellTypes`.
}
\description{
Generate reference data from single-cell expression with cell type label for each cell.
}
\examples{
library(DeconBenchmark)
data(BloodExample) # Load bulk data
print(names(BloodExample)) # c("bulk", "singleCellExpr", "singleCellLabels")
reference <- generateReference(BloodExample$singleCellExpr,
                               BloodExample$singleCellLabels,
                               c("markers", "sigGenes", "signature", "cellTypeExpr"))
print(reference$markers)
print(reference$sigGenes)
print(head(reference$signature))
print(head(reference$cellTypeExpr))

}
