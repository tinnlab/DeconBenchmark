library(rhdf5)

#’ Get supported deconvolution methods
#’
#’ Get supported deconvolution methods
#’
#’ @return a array of supported methods.
#’
#' @examples
#' supportedMethods <- getSupportedMethods()
#’ @export
getSupportedMethods <- function(){
  c("AdRoit", "ARIC", "AutoGeneS", "BayCount", "BayesCCE", "BayICE", "BisqueMarker", "CellDistinguisher",
    "CIBERSORT", "CPM", "DAISM", "debCAM", "Deblender", "DeCompress", "deconf", "DeconICA", "DeconPeaker",
    "DeconRNASeq", "deconvSeq", "DecOT", "DeMixT", "DESeq2", "DSA", "dtangle", "DWLS", "EMeth", "EPIC",
    "FARDEEP", "LinDeconSeq", "Linseed", "MCPcounter", "MethylResolver", "MIXTURE", "MOMF", "MuSic",
    "NITUMID", "PREDE", "quanTIseq", "ReFACTor", "RNA-Sieve", "scaden", "SCDC", "TOAST")
}

#’ Get docker image name
#’
#’ Get docker image name. This method is not to be called by users.
#’
#’ @param method Method name
#’
#’ @return docker image name for the method
.getMethodDockerRepos <- function(method){
  paste0("deconvolution/", tolower(method), ":latest")
}

#’ Get required inputs for a list of methods
#’
#’ Get required inputs for a list of methods
#’
#’ @param methods An array of method names
#’
#’ @return a list with names/keys are methods and values of each key are required parameters of the coressponding method.
#’
#' @examples
#' \dontrun{
#' requiredInputs <- getRequiredInputs(c("AdRoit", "ARIC"))
#' print(requiredInputs$AdRoit)
#' print(requiredInputs$ARIC)
#' }
#’ @importFrom processx run
#’ @export
getMethodsInputs <- function(methods, verbose = T) {
  methodParams <- list()

  for (method in methods) {
    tmpDir <- file.path(tempdir(), paste0(sample(c(LETTERS, letters), 10, TRUE), collapse = ""))

    params <- c("run", "--rm",
                "-v", paste0(tmpDir, ":", "/output"),
                '-e', 'PARAMS_OUTPUT_PATH=/output/params.csv',
                .getMethodDockerRepos(method)
    )

    if (verbose) message("Running docker ", paste(params, collapse = " "))

    tryCatch(
      processx::run("docker",
                    params,
                    error_on_status = FALSE,
                    stdout_callback = function(newout, proc) {
                      if (verbose) message(newout)
                    },
                    stderr_callback = function(newerr, proc) {
                      message(newerr)
                    }),
      error = function(e) {
        list(status = T)
      }
    )

    resultFile <- file.path(tmpDir, "params.csv")
    if (file.exists(resultFile)) {
      methodParam <- read.csv(resultFile, header = T, stringsAsFactors = F)[, 1]
      methodParams[[method]] <- methodParam
    } else {
      methodParams[[method]] <- NULL
      continue()
    }
  }

  methodParams
}

#’ Write arguments to a h5 file
#’
#’ Write deconvolution arguments to a h5 file. This method is not to be called by users.
#’
#’
.writeArgs <- function(h5file,
                       bulk,
                       nCellTypes = NULL,
                       markers = NULL,
                       isMethylation = F,
                       seed = 1,
                       singleCellExpr = NULL,
                       singleCellLabels = NULL,
                       singleCellSubjects = NULL,
                       cellTypeExpr = NULL,
                       sigGenes = NULL,
                       signature = NULL
) {

  writeListOrVector <- function(v, name) {
    h5createGroup(h5file, name)
    h5write(v, h5file, paste0(name, "/values"))
    if (!is.null(names(v))) {
      h5write(names(v), h5file, paste0(name, "/names"))
    }
  }

  writeMatrix <- function(m, name) {
    h5createGroup(h5file, name)
    h5write(m, h5file, paste0(name, "/values"))
    if (!is.null(rownames(m))) {
      h5write(rownames(m), h5file, paste0(name, "/rownames"))
    }
    if (!is.null(colnames(m))) {
      h5write(colnames(m), h5file, paste0(name, "/colnames"))
    }
  }

  unlink(h5file)
  h5createFile(h5file)

  if (!is.null(seed))               writeListOrVector(seed, "seed")
  if (!is.null(bulk))               writeMatrix(bulk, "bulk")
  if (!is.null(nCellTypes))         writeListOrVector(nCellTypes, "nCellTypes")
  if (!is.null(markers))            writeListOrVector(markers, "markers")
  if (!is.null(cellTypeExpr))       writeMatrix(cellTypeExpr, "cellTypeExpr")
  if (!is.null(sigGenes))           writeListOrVector(sigGenes, "sigGenes")
  if (!is.null(signature))          writeMatrix(signature, "signature")
  if (!is.null(singleCellExpr))     writeMatrix(singleCellExpr, "singleCellExpr")
  if (!is.null(singleCellLabels))   writeListOrVector(singleCellLabels, "singleCellLabels")
  if (!is.null(singleCellSubjects)) writeListOrVector(singleCellSubjects, "singleCellSubjects")
  if (!is.null(isMethylation))      writeListOrVector(isMethylation, "isMethylation")
}
