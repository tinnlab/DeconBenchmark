#' @title Get supported deconvolution methods
#'
#' @description Get supported deconvolution methods
#'
#' @return an array of supported methods.
#' @examples
#' supportedMethods <- getSupportedMethods()
#' @export
getSupportedMethods <- function(){
  c("AdRoit", "ARIC", "AutoGeneS", "BayCount", "BayesCCE", "BayICE", "BisqueMarker", "CellDistinguisher",
    "CIBERSORT", "CPM", "DAISM", "debCAM", "Deblender", "DeCompress", "deconf", "DeconICA", "DeconPeaker",
    "DeconRNASeq", "deconvSeq", "DecOT", "DeMixT", "DESeq2", "DSA", "dtangle", "DWLS", "EMeth", "EPIC",
    "FARDEEP", "LinDeconSeq", "Linseed", "MCPcounter", "MethylResolver", "MIXTURE", "MOMF", "MuSic",
    "NITUMID", "PREDE", "quanTIseq", "ReFACTor", "RNA-Sieve", "scaden", "SCDC", "TOAST")
}

#' @title Get docker image name
#'
#' @description Get docker image name. This method is not to be called by users.
#'
#' @param method Method name
#'
#' @return docker image name for the method
.getMethodDockerRepos <- function(method){
  paste0("deconvolution/", tolower(method), ":latest")
}

#' @title Get required inputs for a list of methods
#'
#' @description Get required inputs for a list of methods
#'
#' @param methods An array of method names
#' @param containerEngine A string indicating the container engine, must be `docker` or `singularity`. The default value is `docker`.
#' @param verbose A logical value (TRUE or FALSE) indicating whether to print the output of running processes.
#'
#' @return a list with names/keys are methods and values of each key are required parameters of the coressponding method.
#'
#' @examples
#' \dontrun{
#' requiredInputs <- getRequiredInputs(c("AdRoit", "ARIC"))
#' print(requiredInputs$AdRoit)
#' print(requiredInputs$ARIC)
#' }
#' @importFrom processx run
#' @export
getMethodsInputs <- function(methods, containerEngine = "docker", verbose = T) {
  methodParams <- list()

  for (method in methods) {
    tmpDir <- file.path(tempdir(), paste0(sample(c(LETTERS, letters), 10, TRUE), collapse = ""))
    dir.create(tmpDir, recursive = TRUE, showWarnings = F)

    if (containerEngine == "docker"){
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
    } else if (containerEngine == "singularity"){
      params <- c("run",
                  '--env', paste0('PARAMS_OUTPUT_PATH=', tmpDir, '/params.csv'),
                  paste0("docker://", .getMethodDockerRepos(method))
      )

      if (verbose) message("Running singularity ", paste(params, collapse = " "))

      tryCatch(
        processx::run("singularity",
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
    } else {
      stop("Unknown container engine", containerEngine)
    }

    resultFile <- file.path(tmpDir, "params.csv")
    if (file.exists(resultFile)) {
      methodParam <- read.csv(resultFile, header = T, stringsAsFactors = F)[, 1]
      methodParams[[method]] <- methodParam
    } else {
      methodParams[[method]] <- NULL
    }
    unlink(tmpDir, recursive = TRUE)
  }

  methodParams
}

#' @title Write arguments to a h5 file
#'
#' @description Write deconvolution arguments to a h5 file. This method is not to be called by users.
#' @importFrom rhdf5 h5createGroup h5write h5createFile
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
    rhdf5::h5createGroup(h5file, name)
    rhdf5::h5write(v, h5file, paste0(name, "/values"))
    if (!is.null(names(v))) {
      rhdf5::h5write(names(v), h5file, paste0(name, "/names"))
    }
  }

  writeMatrix <- function(m, name) {
    rhdf5::h5createGroup(h5file, name)
    rhdf5::h5write(m, h5file, paste0(name, "/values"))
    if (!is.null(rownames(m))) {
      rhdf5::h5write(rownames(m), h5file, paste0(name, "/rownames"))
    }
    if (!is.null(colnames(m))) {
      rhdf5::h5write(colnames(m), h5file, paste0(name, "/colnames"))
    }
  }

  unlink(h5file)
  rhdf5::h5createFile(h5file)

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
