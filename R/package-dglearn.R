#' Learning with Diffusion Geometry
#' 
#' @docType package
#' @noRd
#' @import Rdpack
#' @import randomForest
#' @importFrom igraph distances graph_from_adjacency_matrix
#' @importFrom stats rnorm dist as.dist
#' @importFrom Rcpp evalCpp
#' @useDynLib dglearn
NULL
# pack <- "repsim"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))