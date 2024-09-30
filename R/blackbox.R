#' Blackbox Scaling
#'
#' \code{blackbox} is a function that takes a matrix of survey data in which
#' individuals place themselves on continuous scales across multiple issues,
#' and locates those citizens in a spatial model of voting. Mathematically,
#' this function generalizes the singular value of a matrix to cases in which
#' there is missing data in the matrix. Scales generated using perceptual data
#' (i.e. scales of legislator locations using liberal-conservative rankings by
#' survey respondents) should instead use the \code{blackbox_transpose}
#' function in this package instead.
#'
#'
#' @param data matrix of numeric values containing the issue scale data.
#' Respondents should be organized on rows, and stimuli on columns. It is
#' helpful, though not necessary, to include row names and column names.
#' @param missing vector or matrix of numeric values, sets the missing values
#' for the data.  NA values are always treated as missing regardless of what is
#' set here.  Observations with missing data are discarded before analysis.  If
#' input is a vector, then the vector is assumed to contain the missing value
#' codes for all the data.  If the input is a matrix, it must be of dimension p
#' x q, where p is the maximum number of missing values and q is the number of
#' columns in the data.  Each column of the inputted matrix then specifies the
#' missing data values for the respective variables in data.  If null
#' (default), no missing values are in the data other than the standard NA
#' value.
#' @param verbose logical, indicates whether \code{aldmck} should print out
#' detailed output when scaling the data.
#' @param dims integer, specifies the number of dimensions to be estimated.
#' @param minscale integer, specifies the minimum number of responses a
#' respondent needs needs to provide to be used in the scaling.
#' @return An object of class \code{blackbox}.
#'
#' \item{stimuli}{ vector of data frames of length dims. Each data frame
#' presents results for estimates from that dimension (i.e. `x$stimuli[[2]]`
#' presents results for dimension 2).  Each row contains data on a separate
#' stimulus, and each data frame includes the following variables: \itemize{
#' \item\code{N}Number of respondents who provided a response to this stimulus.
#' \item\code{c}Stimulus intercept.  \item\code{w1}Estimate of the stimulus
#' weight on the first dimension. If viewing the results for a higher
#' dimension, higher dimension results will appear as w2, w3, etc.
#' \item\code{R2}The percent variance explained for the stimulus. This
#' increases as more dimensions are estimated.  } }
#'
#' \item{individuals}{ vector of data frames of length dims. Each data frame
#' presents results for estimates from that dimension (i.e. `x$stimuli[[2]]`
#' presents results for dimension 2).  Individuals that are discarded from
#' analysis due to the minscale constraint are NA'd out.  Each row contains
#' data on a separate stimulus, and each data frame includes the following
#' variables: \itemize{ \item\code{c1}Estimate of the individual intercept on
#' the first dimension. If viewing the results for a higher dimension, higher
#' dimension results will appear as c2, c3, etc.  } } \item{fits}{ A data frame
#' of fit results, with elements listed as follows:} \itemize{
#' \item\code{SSE}Sum of squared errors.  \item\code{SSE.explained}Explained
#' sum of squared error.  \item\code{percent}Percentage of total variance
#' explained.  \item\code{SE}Standard error of the estimate, with formula
#' provided on pg. 973 of the article cited below.
#' \item\code{singular}Singluar value for the dimension.  } \item{Nrow}{ Number
#' of rows/stimuli.} \item{Ncol}{ Number of columns used in estimation. This
#' may differ from the data set due to columns discarded due to the minscale
#' constraint.} \item{Ndata}{ Total number of data entries.} \item{Nmiss}{
#' Number of missing entries.} \item{SS_mean}{ Sum of squares grand mean.}
#' \item{dims}{ Number of dimensions estimated.}
#' @seealso '\link{Issues1980}'.
#' @references
#'
#' Keith Poole, Jeffrey Lewis, Howard Rosenthal, James Lo, Royce Carroll (2016)
#' ``Recovering a Basic Space from Issue Scales in R.'' Journal of Statistical
#' Software. 69(7), 1--21. doi:10.18637/jss.v069.i07
#'
#' Keith T. Poole (1998) ``Recovering a Basic Space From a Set of Issue
#' Scales.'' American Journal of Political Science. 42(3), 954-993.
#' @keywords multivariate
#' @export
#' @examples
#'
#'
#' ### Loads issue scales from the 1980 NES.
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#' summary(Issues1980_bb)
#'
#'
blackbox <- function(data,missing=NULL,verbose=FALSE,dims=1,minscale){

  ### Error check each argument ###
  if(inherits(data, "data.frame"))  data <- as.matrix(data)
  if(!inherits(data, "matrix"))  stop("Data is not a matrix, or data frame.")

  if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric() or mode().")
  if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
  if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")
  if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")
  if(minscale<1) stop("Argument 'minscale' must be positive.")
  if(dims<1) stop("Argument 'dims' must be positive.")

  ### Format the data input ###
  N <- nrow(data)
  NQ <- ncol(data)

  if(is.vector(missing))	data[data %in% missing] <- NA
  if(is.matrix(missing))	for(i in 1:ncol(data))	data[data[,i] %in% missing[,i],i] <- NA
  missval <- max(data,na.rm=TRUE) + 1	#code for missing data
  rawdata <- as.numeric(t(data))
  rawdata[is.na(rawdata)] <- missval

  ##Longer output
  stimnames <- colnames(data)
  if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
  if(verbose){

    deleted <- sum(is.na(apply(data,1,sum)))
    cat("\n\n\tBeginning Blackbox Scaling...")

    #cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
    cat(NQ, "stimuli have been provided.")

  }

  res <- blackboxf_wrapper(
                  as.integer(N),			# NRESPONDENTS
                  as.integer(NQ),			# NISSUES
                  as.integer(dims),		# NDIMENSIONS, check later about mods
                  as.integer(1),  		# NMISSING
                  as.numeric(rep(missval,NQ)),	# KMISS
                  as.integer(minscale),		# MINSCALE
                  as.integer(rep(1,N)),		# MID
                  as.numeric(rawdata),		# KISSUE
                  fits = double(7*dims), 		# FITS
                  psimatrix = double(N*((dims*(dims+1))/2)+2*N*dims),		# PSIMATRIX
                  wmatrix = double((NQ)*((dims*(dims+1))/2)+2*(NQ)*dims),		# WMATRIX
                  lresp = integer(N+NQ),						# LRESPONDENTS
                  lmark = integer(N),		# LMARK
                  fits2 = double(6),		# FITS2
                  exitstatus = integer(1))	# EXITSTATUS

  if (res$exitstatus != 1) stop("\n\n\t====== Blackbox did not execute properly ======\n\n")

  stimuli <- vector("list", dims)
  start <- 1
  end <- 3*NQ
  for(i in 1:dims){
    stimuli[[i]] <- as.data.frame(matrix(round(res$wmatrix[start:end],digits=3),nrow=NQ,ncol=i+2,byrow=T))
    colnames(stimuli[[i]]) <- c("c",paste("w",1:i,sep=""),"R2")
    rownames(stimuli[[i]]) <- stimnames
    stimuli[[i]] <- cbind(N=res$lresp[1:NQ],stimuli[[i]])
    start <- end + 1
    end <- start + (i+3)*NQ - 1
  }
  individuals <- vector("list", dims)
  start <- 1
  end <- N
  for(i in 1:dims){
    individuals[[i]] <- as.data.frame(matrix(round(res$psimatrix[start:end],digits=3),nrow=N,ncol=i,byrow=T))
    dumpthese <- (rowSums(individuals[[i]]==0)==i)
    individuals[[i]][dumpthese,] <- NA
    colnames(individuals[[i]]) <- c(paste("c",1:i,sep=""))
    rownames(individuals[[i]]) <- rownames(data)
    start <- end + 1
    end <- start + (i+1)*N - 1
  }

  fits <- matrix(res$fits,nrow=dims,ncol=7,byrow=T)
  fits <- as.data.frame(fits[,c(1:3,6:7),drop=FALSE])
  colnames(fits) <- c("SSE","SSE.explained","percent","SE","singular")
  rownames(fits) <- paste("Dimension",1:dims)

  result <- list(stimuli = stimuli, individuals=individuals, fits=fits,
                 Nrow = res$fits2[1],
                 Ncol = res$fits2[2],
                 Ndata = res$fits2[3],
                 Nmiss = res$fits2[4],
                 SS_mean = res$fits2[6],
                 dims=dims)

  class(result) <- c("blackbox")
  if (verbose) cat("\n\n\tBlackbox estimation completed successfully.\n\n")
  result

}

#' Blackbox Summary
#'
#' \code{summary.blackbox} reads an \code{blackbox} object and prints a
#' summary.
#'
#'
#' @param object a \code{blackbox} output object.
#' @param ...  further arguments to \code{print}.
#' @return A summary of a \code{blackbox} object. For each dimension, reports
#' all stimuli with coordinates, individuals used for scaling, and fit. Also
#' summarizes number of rows, columns, total data entries, number of missing
#' entries, percent missing data, and sum of squares.
#' @seealso '\link{blackbox}', '\link{Issues1980}'
#' @keywords multivariate
#' @export
#' @method summary blackbox
#' @examples
#'
#' ### Loads issue scales from the 1980 NES.
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#' summary(Issues1980_bb)
#'
summary.blackbox <- function(object, ...){

  x <- object
  if(!inherits(x, "blackbox"))  stop("Input is not of class blackbox.")

  cat("\n\nSUMMARY OF BLACKBOX OBJECT")
  cat("\n----------------------------------")
  for(i in 1:x$dims){
    cat("\n")
    print(x$stimuli[[i]], ...)
  }

  cat("\n\tDimensions Estimated:", x$dims)
  cat("\n\tNumber of Rows:", x$Nrow)
  cat("\n\tNumber of Columns:", x$Ncol)
  cat("\n\tTotal Number of Data Entries:", x$Ndata)
  cat("\n\tNumber of Missing Entries:", x$Nmiss)
  cat("\n\tPercent Missing Data: ", round(100*x$Nmiss/(x$Ndata + x$Nmiss),2),"%",sep="")
  cat("\n\tSum of Squares (Grand Mean):", x$SS_mean)
  cat("\n\n")
}

#' Blackbox Coordinate Distribution Plot
#'
#' \code{plot.blackbox} reads an \code{blackbox} object and plots a histogram
#' of the estimated intercepts.
#'
#'
#' @param x an \code{blackbox} output object.
#' @param ...  other arguments to \code{hist}.
#' @return
#'
#' A histogram of the estimated intercepts.
#' @seealso \link{Issues1980}.
#' @keywords multivariate
#' @export
#' @method plot blackbox
#' @examples
#'
#'
#'
#' ### Loads issue scales from the 1980 NES.
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#'
#' plot(Issues1980_bb)
#'
#'
plot.blackbox <- function(x, ...){

  if(!inherits(x, "blackbox"))  stop("Input is not of class blackbox.")

  hist(x$individuals[[x$dims]]$c1, breaks = seq(-1, 1, 0.1),
       main="Distribution of Blackbox Intercepts",
       xlab="Blackbox Intercept", ...)

}


#' Predict method of blackbox objects
#'
#' \code{predict.blackbox} reads an \code{blackbox} object and uses the
#' estimates to generate a matrix of predicted values.
#'
#'
#' @param object A \code{blackbox} output object.
#' @param dims Number of dimensions used in prediction. Must be equal to or
#' less than number of dimensions used in estimation.
#' @param ...  Ignored.
#' @return A matrix of predicted values generated from the parameters estimated
#' from a \code{blackbox} object.
#' @seealso '\link{blackbox}', '\link{Issues1980}'
#' @keywords multivariate
#' @export
#' @method predict blackbox
#' @examples
#'
#' ## Estimate blackbox object from example and call predict function
#' data(Issues1980)
#' Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
#' Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
#'
#' ### This command conducts estimates, which we instead load using data()
#' # Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
#' data(Issues1980_bb)
#' prediction <- predict(Issues1980_bb,dims=3)
#'
#' ## Examine predicted vs. observed values for first 10 respondents
#' ## Note that 4th and 6th respondents are NA because of missing data
#' Issues1980[1:10,]
#' prediction[1:10,]
#'
#' ## Check correlation across all predicted vs. observed, excluding missing values
#' prediction[which(Issues1980 %in% c(0,8,9))] <- NA
#' cor(as.numeric(prediction), as.numeric(Issues1980), use="pairwise.complete")
#'
predict.blackbox <- function(object, dims=1, ...){

  ## Error catch for dims
  if(!inherits(object, "blackbox"))  stop("Input is not of class blackbox.")
  if(dims < 1)  stop("dims must be great than 1")
  if(dims > object$dims)  stop(paste("dims must be equal or less than the number of estimate dimensions, which is", object$dims))

  ## Extract object output
  W.hat<- as.matrix(object$stimuli[[dims]][,paste("w", 1:dims, sep="")])
  Psi.hat <- as.matrix(object$individuals[[dims]][,paste("c", 1:dims, sep="")])
  Jn <- rep(1, nrow(Psi.hat))
  c <- object$stimuli[[dims]][,"c"]

  ## In blackbox, Keith already postmultiplied by sqrt(singular), so no need to do it here
  ## This is NOT the case in blackbox_transpose, which needs to be multiplied by sqrt(singular)
  X.hat <- Psi.hat %*% t(W.hat) + Jn %o% c

  colnames(X.hat) <- rownames(object$stimuli[[dims]])
  rownames(X.hat) <- rownames(object$individuals[[dims]])
  return(X.hat)
}





