#' @title Calculate markers from single cell data
#'
#' @description Calculate markers from single cell data. This method is not to be called by users.
#' @param singleCellExpr A matrix of single cell expression data.
#' @param singleCellLabels A vector of cell type labels.
#' @param log2Threshold The log2 threshold for marker selection.
#'
#' @return A list with names/keys are cell type labels and values are the corresponding markers.
#'
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
.generateReference.markers <- function(singleCellExpr, singleCellLabels, log2Threshold = 1) {
  singleCellLabels <- as.character(singleCellLabels)
  colnames(singleCellExpr) <- singleCellLabels

  #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
  keeps <- sapply(unique(singleCellLabels), function(cellType) {
    hits <- singleCellLabels == cellType
    rowSums(singleCellExpr[, hits, drop = FALSE] != 0) >= ceiling(0.3 * sum(hits))
  })

  #normalization
  singleCellExpr <- singleCellExpr[rowSums(keeps) > 0,]
  singleCellExpr <- edgeR::DGEList(singleCellExpr)
  singleCellExpr <- edgeR::calcNormFactors(singleCellExpr, method = "TMM")

  # INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account
  #[compare one group with average expression of all other groups]
  annotation <- factor(singleCellLabels)
  design <- model.matrix(~0 + annotation)
  colnames(design) <- sapply(strsplit(colnames(design), "annotation"), function(x) x[2])
  constrastMatrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
  colnames(constrastMatrix) <- colnames(design)
  diag(constrastMatrix) <- (ncol(design) - 1) / ncol(design)

  v <- limma::voom(singleCellExpr, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, constrastMatrix)
  fit <- limma::eBayes(fit, trend = TRUE)

  topTableResults <- limma::topTable(fit, coef = seq_len(ncol(constrastMatrix)), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2Threshold)
  topTableResults <- topTableResults[, 1:(ncol(topTableResults) - 4)]

  ERCCGenes <- grep("ERCC-", topTableResults$gene)
  if (length(ERCCGenes) > 0) {
    topTableResults <- topTableResults[-ERCCGenes,]
  }

  markers <- apply(topTableResults, 1, function(x) {
    temp <- sort(x)
    ((temp[ncol(topTableResults)] - temp[ncol(topTableResults) - 1]) >= log2Threshold) |
      (abs(temp[1] - temp[2]) >= log2Threshold)
  })

  topTableResults <- topTableResults[markers,]

  markers <- cbind.data.frame(
    rownames(topTableResults),
    t(apply(topTableResults, 1, function(x) {
      temp <- max(x)
      if (temp < log2Threshold) {
        temp <- c(min(x), colnames(topTableResults)[which.min(x)])
      } else {
        temp <- c(max(x), colnames(topTableResults)[which.max(x)])
      }
      temp
    }))
  )

  colnames(markers) <- c("gene", "log2FC", "cellType")
  genes <- as.character(markers$gene)
  cellTypes <- as.character(markers$cellType)

  markers <- lapply(unique(cellTypes), function(cellType) genes[cellTypes == cellType])
  names(markers) <- unique(cellTypes)

  markers
}

#' @title Generate reference data
#'
#' @description Generate reference data from single-cell expression with cell type label for each cell.
#'
#' @param singleCellExpr Single-cell expression matrix in non-log scale. Rows are genes and columns are cells.
#' @param singleCellLabels An array contains cell type label for each cell. The order of labels must be the same as the order of cells in in singleCellExpr.
#' @param types An array of different type of reference to generate.
#' @param log2Threshold The threshold for log2FC for filtering markers.
#'
#' @return A list with names/keys are the same as in types.
#' If type is `markers`, the value is a list with names/keys are unique cell types and values for each key are an array of markers for the corresponding cell type.
#' If type is `sigGenes`, the value is the unique concatenation of all `markers`.
#' If type is `signature`, the value is an matrix of `sigGenes` x `cellTypes`.
#' If type is `cellTypeExpr`,  the value is an matrix of `all genes` x `cellTypes`.
#'
#' @examples
#' library(DeconBenchmark)
#' data(BloodExample) # Load bulk data
#' print(names(BloodExample)) # c("bulk", "singleCellExpr", "singleCellLabels")
#' reference <- generateReference(BloodExample$singleCellExpr,
#'                                BloodExample$singleCellLabels,
#'                                c("markers", "sigGenes", "signature", "cellTypeExpr"))
#' print(reference$markers)
#' print(reference$sigGenes)
#' print(head(reference$signature))
#' print(head(reference$cellTypeExpr))
#'
#' @export
generateReference <- function(singleCellExpr, singleCellLabels, types = c("markers", "sigGenes", "signature", "cellTypeExpr"), log2Threshold = 1) {
  reference <- list()

  if (any(c("markers", "sigGenes", "signature") %in% types)) {
    reference$markers <- .generateReference.markers(singleCellExpr, singleCellLabels, log2Threshold)
  }

  if (any(c("sigGenes", "signature") %in% types)) {
    reference$sigGenes <- unique(unlist(reference$markers))
  }

  if (any(c("signature", "cellTypeExpr") %in% types)) {
    cellTypeExpr <- matrix(ncol = length(unique(singleCellLabels)), nrow = nrow(singleCellExpr))
    colnames(cellTypeExpr) <- unique(singleCellLabels)
    rownames(cellTypeExpr) <- rownames(singleCellExpr)

    for (cellType in unique(singleCellLabels)) {
      tmp <- rowSums(singleCellExpr[, singleCellLabels == cellType])
      tmp <- tmp / sum(tmp) * 1e6
      cellTypeExpr[, cellType] <- tmp
    }

    reference$cellTypeExpr <- cellTypeExpr
  }

  if ("signature" %in% types) {
    reference$signature <- reference$cellTypeExpr[reference$sigGenes, ]
  }

  reference[types]
}
