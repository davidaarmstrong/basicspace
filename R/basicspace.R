#' basicspace: Recovering a Basic Space from Issue Scales.
#'
#' Conducts Aldrich-McKelvey and Blackbox Scaling (Poole et al 2016)
#' <doi:10.18637/jss.v069.i07> to recover latent dimensions of judgment.
#'
#' @section basicspace functions:
#' The main functions are `aldmck()` and `blackbox()` and `blackbox_transpose()`.
#'
#' @import tools
#' @importRcpp
#' @importFrom grDevices palette
#' @importFrom graphics abline arrows hist mtext par plot points segments text
#' @importFrom stats cor density ecdf quantile
#' @importFrom utils flush.console
#'
## usethis namespace: start
#' @useDynLib basicspace, .registration = TRUE
## usethis namespace: end
#'
#' @name basicspace
"_PACKAGE"

