#' Stimulus extraction function
#'
#' \code{stimuli} is a convenience function to extract the stimulus parameters
#' from an \code{aldmck}, \code{blackbox}, or \code{blackbt} object.
#'
#'
#' @param object an \code{aldmck}, \code{blackbox}, or \code{blackbt} output
#' object.
#' @return The stimuli of the estimated output, which can also be recovered as
#' \code{object$stimuli}.  Please refer to the documentation of \code{aldmck},
#' \code{blackbox}, or \code{blackbox_transpose} for specifics.
#' @seealso '\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
#' @keywords multivariate
#' @export
#' @examples
#'
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#' stimuli(Issues1980_bb)
#'
stimuli <- function(object){

  if(inherits(object, "blackbox") | inherits(object, "blackbt")) return(object$stimuli)
  if(inherits(object, "aldmck") )  return(object$stimuli)
  stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")

}

#' Extraction function for scaled individuals
#'
#' \code{individuals} is a convenience function to extract the
#' individual/respondent parameters from an \code{aldmck}, \code{blackbox}, or
#' \code{blackbt} object.
#'
#'
#' @param object an \code{aldmck}, \code{blackbox}, or \code{blackbt} output
#' object.
#' @return The individual parameters of the estimated output, which can also be
#' recovered as \code{object$individuals} (for \code{blackbox} or
#' \code{blackbt} objects) or \code{object$respondents} (for \code{aldmck}
#' objects).  Please refer to the documentation of \code{aldmck},
#' \code{blackbox}, or \code{blackbox_transpose} for specifics.
#' @seealso '\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
#' @keywords multivariate
#' @export
#' @examples
#'
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#' individuals(Issues1980_bb)
#'
individuals <- function(object){

  if(inherits(object, "blackbox") | inherits(object, "blackbt")) return(object$individuals)
  if(inherits(object, "aldmck") )  return(object$respondents)
  stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")

}



#' Extraction function for fit of scaling model
#'
#' \code{fit} is a convenience function to extract the model fit statistics
#' from an \code{aldmck}, \code{blackbox}, or \code{blackbt} object.
#'
#'
#' @param object an \code{aldmck}, \code{blackbox}, or \code{blackbt} output
#' object.
#' @return The model fit statistics of the estimated output, which can also be
#' recovered as \code{object$fits} (for \code{blackbox} or \code{blackbt}
#' objects) or \code{object$AMfit} (for \code{aldmck} objects).  Please refer
#' to the documentation of \code{aldmck}, \code{blackbox}, or
#' \code{blackbox_transpose} for specifics.
#' @seealso '\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
#' @keywords multivariate
#' @export
#' @examples
#'
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#'
#'
#' fit(Issues1980_bb)
#'
#'
fit <- function(object){

  if(inherits(object, "blackbox") | inherits(object, "blackbt")) return(object$fits)
  if(inherits(object, "aldmck") )  return(object$AMfit)
  stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")
}

#' @method dim blackbox
#' @export
dim.blackbox <- function(x){
  return(c(x$Nrow, x$Ncol))
}

#' @method dim blackbt
#' @export
dim.blackbt <- function(x){
  return(c(x$Nrow, x$Ncol))
}

#' @method dim aldmck
#' @export
dim.aldmck <- function(x){
  return(c(x$N, length(x$stimuli)))
}

#' @method ncol blackbox
#' @export
ncol.blackbox <- function(x){
  return(x$Ncol)
}

#' @method nrow blackbox
#' @export
nrow.blackbox <- function(x){
  return(x$Nrow)
}

#' @method ncol blackbt
#' @export
ncol.blackbt <- function(x){
  return(x$Ncol)
}

#' @method nrow blackbt
#' @export
nrow.blackbt <- function(x){
  return(x$Nrow)
}

#' @method nrow aldmck
#' @export
nrow.aldmck <- function(x){
  return(x$N)
}

#' @method ncol aldmck
#' @export
ncol.aldmck <- function(x){
  return(length(x$stimuli))
}

#' Coordinate Cumulative Distribution Plot
#'
#' Plots the cumulative distribution of the respondents and stimuli.
#'
#'
#' @param x an \code{aldmck} output object.
#' @param align integer, the x-axis location that stimuli names should be
#' aligned to If set to NULL, it will attempt to guess a location.
#' @param xlim vector of length 2, fed to the \code{plot} function as the
#' \code{xlim} argument, which sets the minimum and maximum range of the
#' x-axis.
#' @param ...  other arguments to \code{plot}.
#' @return A plot of the empirical cumulative distribution of the respondent
#' ideal points, along with the locations of the stimuli. If no self-placements
#' were specified during estimation, no graphical plots will appear.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}',
#' '\link{plot.aldmck}'.
#' @keywords multivariate
#' @export
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#'
#' summary(result)
#' par(ask=TRUE)
#' plotcdf(result)
#'
plotcdf <- function(x, ..., align=NULL, xlim=c(-2,2)){
  UseMethod("plotcdf")
}
