#’ Run deconvolution
#’
#’ Run deconvolution for a list of methods.
#’ See getSupportedMethods() for the list of supported methods.
#  Each method is run in a separate docker container.
#’ Each method has its own required inputs.
#  See getMethodsInputs() for the input of each method.
#’
#’ @param methods An array of method names. See getSupportedMethods() for the list of supported methods.
#’ @param bulk A matrix of bulk data. Rows are genes and columns are bulk samples.
#’ @param nCellTypes The number of cell types.
#’ @param markers A list with name/key are cell types and values of each key are markers for the coressponding cell type.
#’ @param isMethylation A logical value (TRUE or FALSE) indicating whether the input data is methylation data.
#’ @param singleCellExpr A matrix of single cell data. Rows are genes and columns are single cell samples.
#’ @param singleCellLabels A array of single cell labels. The order of the labels must be the same as the order of the single cell samples.
#’ @param singleCellSubjects A array of single cell subjects. The order of the subjects must be the same as the order of the single cell samples.
#’ @param cellTypeExpr A matrix of cell type data. Rows are genes and columns are cell types. This matrix normally contains the same number of genes as bulk expression matrix (not to be confused with signature matrix).
#’ @param sigGenes A list of significant genes.
#’ @param signature A matrix of signature data. Rows are genes and columns are cell types. This matrix normally contains only marker genes (not to be confused with cell type expresison matrix).
#’ @param seed An integer value for the seed of the random number generator.
#’ @param matlabLicenseFile A path to a matlab license file. This file is required for methods that use matlab. The hostid of this license is the same as the hostid of the host machine. The user of this license must be `root`. Visit Matlab license center (https://www.mathworks.com/licensecenter/licenses) to obtain a license file.
#’ @param timeout A timeout in seconds for the docker container. The default value is 12 hours.
#’ @param dockerArgs A list of extra arguments for the docker container. The default values allow docker container to use maximum 8 cpus, reserve 4GB of memory and allow it to use up to 32GB of memory.
#’ @param verbose A logical value (TRUE or FALSE) indicating whether to print the output of running processes.
#’
#’ @return A list with names/keys are methods and values of each key is the deconvolution result for the corresponding methods. The deconvolution result for each method is a list with two elements: `$P` and `$S`. The `$P` is a matrix cell type proportions (samples x cell type) and the `$S` is a matrix of cell type signatures (genes x cell type).
#’
#’ @examples
#’ \dontrun{
#’ library(DeconUtils)
#’ allSupportMethods <- getSupportedMethods() # Get the list of supported methods
#’ print(allSupportMethods)
#’ methodsToRun <- c("AdRoit", "ARIC") # Select methods to run (must be in the list of supported methods)
#’ requiredInputs <- getMethodsInputs(supportMethods) # Get the required inputs for each method
#’ print(requiredInputs) # list (AdRoif = c("bulk", "singleCellExpr", "singleCellLabels"), ARIC = c("bulk", "cellTypeExpr", "isMethylation"))
#’ data(BloodExample) # Load example data
#’ print(names(BloodExample)) # c("bulk", "singleCellExpr", "singleCellLabels")
#’ bulk <- BloodExample$bulk
#’ singleCellExpr <- BloodExample$singleCellExpr
#’ singleCellLabels <- BloodExample$singleCellLabels
#’ # Run AdRoit only
#’ deconvolutionResult <- deconvolution(methods = "AdRoit", bulk = bulk, singleCellExpr = singleCellExpr, singleCellLabels = singleCellLabels)
#’ proportion <- deconvolutionResult$AdRoit$P
#’ print(proportion)
#’ # Run ARIC only
#’ reference <- generateReference(singleCellExpr, singleCellLabels, type="cellTypeExpr") # Generate reference
#’ deconvolutionResult <- deconvolution(methods = "ARIC", bulk = bulk,  cellTypeExpr = reference$cellTypeExpr, isMethylation = F)
#’ proportion <- deconvolutionResult$ARIC$P
#’ print(proportion)
#’ # Run both AdRoit and ARIC
#’ reference <- generateReference(singleCellExpr, singleCellLabels, type="cellTypeExpr") # Generate reference
#’ deconvolutionResults <- deconvolution(methodsToRun, bulk = bulk, singleCellExpr = singleCellExpr, singleCellLabels = singleCellLabels, cellTypeExpr = reference$cellTypeExpr, isMethylation = F) # Run deconvolution
#’ proportionAdRoit <- deconvolutionResults$AdRoit$P
#’ proportionARIC <- deconvolutionResults$AdRoit$P
#’ print(proportionAdRoit)
#’ print(proportionARIC)
#’ }
#’ @import processx
#’ @export
runDeconvolution <- function(methods,
                             bulk,
                             nCellTypes = NULL,
                             markers = NULL,
                             isMethylation = F,
                             singleCellExpr = NULL,
                             singleCellLabels = NULL,
                             singleCellSubjects = NULL,
                             cellTypeExpr = NULL,
                             sigGenes = NULL,
                             signature = NULL,
                             seed = 1,
                             matlabLicenseFile = NULL,
                             timeout = 60*60*12,
                             dockerArgs = c(
                               '--cpus=8.0',
                               '-m=32G',
                               '--memory-reservation=4G'
                             ),
                             verbose = T) {

  tmpDir <- file.path(tempdir(), paste0(sample(c(LETTERS, letters), 10, TRUE), collapse = ""))
  dir.create(tmpDir, recursive = TRUE, showWarnings = F)
  tmpH5File <- file.path(tmpDir, "args.h5")

  .writeArgs(h5file = tmpH5File,
             bulk = bulk,
             nCellTypes = nCellTypes,
             markers = markers,
             isMethylation = isMethylation,
             singleCellExpr = singleCellExpr,
             singleCellLabels = singleCellLabels,
             singleCellSubjects = singleCellSubjects,
             cellTypeExpr = cellTypeExpr,
             sigGenes = sigGenes,
             signature = signature,
             seed = seed)

  .message <- function(...) {
    if (verbose) message(paste(...))
  }

  allResults <- list()

  for (method in methods) {

    params <- c("run", "--rm", dockerArgs,
                "-v", paste0(tmpH5File, ":", "/args.h5"),
                "-v", paste0(tmpDir, ":", "/output"),
                '-e', paste0('OUTPUT_PATH=/output/', method,'-results.h5')
    )

    if (!is.null(matlabLicenseFile)) {
      params <- c(params, "-v", paste0(matlabLicenseFile, ':', '/licenses/license.lic'), '--network=host')
    }
    params <- c(params, .getMethodDockerRepos(method))

    .message("Running docker ", paste(params, collapse = " "))

    res <- tryCatch(
      processx::run("docker",
                    params,
                    timeout = timeout,
                    error_on_status = FALSE,
                    stdout_callback = function(newout, proc) {
                      .message(newout)
                    },
                    stderr_callback = function(newerr, proc) {
                      message(newerr)
                    }),
      error = function(e) {
        list(status = T)
      }
    )

    if (res$status) {
      allResults[[method]] <- NULL
      next()
    }

    resFile <- file.path(tmpDir, paste0(method, "-results.h5"))
    if (length(resFile) == 0) {
      allResults[[method]] <- NULL
      next()
    }

    P <- tryCatch(
      h5read(resFile, "P"),
      error = function(e) {
        NULL
      }
    )

    S <- tryCatch(
      h5read(resFile, "S"),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(P)) {
      .P <- P$values
      if (!is.null(P$rownames)) {
        rownames(.P) <- P$rownames
      }
      if (!is.null(P$colnames)) {
        colnames(.P) <- P$colnames
      }
      P <- .P
    }

    if (!is.null(S)) {
      .S <- S$values
      if (!is.null(S$rownames)) {
        rownames(.S) <- S$rownames
      }
      if (!is.null(S$colnames)) {
        colnames(.S) <- S$colnames
      }
      S <- .S
    }

    allResults[[method]] <- list(P = P, S = S)
  }

  unlink(tmpDir, recursive = TRUE)

  allResults
}