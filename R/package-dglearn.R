#' Learning with Diffusion Geometry
#' 
#' @docType package
#' @noRd
#' @import randomForest
#' @import maotai 
#' @importFrom utils getFromNamespace
#' @importFrom Rdpack reprompt
#' @importFrom RSpectra eigs
#' @importFrom igraph distances graph_from_adjacency_matrix
#' @importFrom stats rnorm dist as.dist median
#' @importFrom Rcpp evalCpp
#' @useDynLib dglearn
NULL
# pack <- "repsim"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))