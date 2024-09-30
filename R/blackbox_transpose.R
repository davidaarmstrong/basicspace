#' Blackbox transpose Scaling
#'
#' \code{blackbox_transpose} is a function that takes a matrix of perceptual
#' data, such as liberal-conservative rankings of various stimuli, and recovers
#' the true location of those stimuli in a spatial model. It differs from
#' procedures such as \code{wnominate}, which instead use preference data to
#' estimate candidate and citizen positions. The procedure here generalizes the
#' technique developed by John Aldrich and Richard McKelvey in 1977, which is
#' also included in this package as the \code{aldmck} function.
#'
#'
#' @param data matrix of numeric values, containing the perceptual data.
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
#' @return An object of class \code{blackbt}.
#'
#' \item{stimuli}{ vector of data frames of length dims. Each data frame
#' presents results for estimates from that dimension (i.e. `x$stimuli[[2]]`
#' presents results for dimension 2).  Each row contains data on a separate
#' stimulus, and each data frame includes the following variables: \itemize{
#' \item\code{N}Number of respondents who ranked this stimulus.
#' \item\code{coord1D}Location of the stimulus in the first dimension. If
#' viewing the results for a higher dimension, higher dimension results will
#' appear as coord2D, coord3D, etc.  \item\code{R2}The percent variance
#' explained for the stimulus. This increases as more dimensions are estimated.
#' } }
#'
#' \item{individuals}{ vector of data frames of length dims. Each data frame
#' presents results for estimates from that dimension (i.e. `x$stimuli[[2]]`
#' presents results for dimension 2).  Individuals that are discarded from
#' analysis due to the minscale constraint are NA'd out.  Each row contains
#' data on a separate stimulus, and each data frame includes the following
#' variables: \itemize{ \item\code{c}Estimate of the individual intercept.
#' \item\code{w1}Estimate of the individual slope. If viewing the results for a
#' higher dimension, higher dimension results will appear as w2, w3, etc.
#' \item\code{R2}The percent variance explained for the respondent. This
#' increases as more dimensions are estimated.  } } \item{fits}{ A data frame
#' of fit results, with elements listed as follows:} \itemize{
#' \item\code{SSE}Sum of squared errors.  \item\code{SSE.explained}Explained
#' sum of squared error.  \item\code{percent}Percentage of total variance
#' explained.  \item\code{SE}Standard error of the estimate, with formula
#' provided in the article cited below.  \item\code{singular}Singluar value for
#' the dimension.  } \item{Nrow}{ Number of rows/stimuli.} \item{Ncol}{ Number
#' of columns used in estimation. This may differ from the data set due to
#' columns discarded due to the minscale constraint.} \item{Ndata}{ Total
#' number of data entries.} \item{Nmiss}{ Number of missing entries.}
#' \item{SS_mean}{ Sum of squares grand mean.} \item{dims}{ Number of
#' dimensions estimated.}
#' @seealso '\link{plotcdf}', '\link{LC1980}', '\link{LC1980_bbt}'.
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
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' LCdat=LC1980[,-1]	#Dump the column of self-placements
#'
#'
#' ### This command conducts estimates, which we instead load using data()
#'
#' #LC1980_bbt <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
#' data(LC1980_bbt)
#' plot(LC1980_bbt)
#'
#' par(ask=TRUE)
#' plotcdf(LC1980_bbt)
#' summary(LC1980_bbt)
#'
blackbox_transpose <- function(data,missing=NULL,verbose=FALSE,dims=1,minscale){

  ### Error check each argument ###
  if(inherits(data, "data.frame"))  data <- as.matrix(data)
  if(!inherits(data, "matrix"))  stop("Data is not a matrix, or data frame.")

  if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric().")
  if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
  if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")
  if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")
  if(minscale<1) stop("Argument 'minscale' must be positive.")
  if(dims<1) stop("Argument 'dims' must be positive.")
  if(nrow(data)>1500) stop("There are more than N = 1500 respondents in the data.")

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
    cat("\n\n\tBeginning Blackbox Transpose Scaling...")

    #cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
    cat(NQ, "stimuli have been provided.")

  }

  res <- blackboxt_wrapper(
                  as.integer(N),			# NRESPONDENTS
                  as.integer(NQ),			# NISSUES
                  as.integer(dims),		# NDIMENSIONS, check later about mods
                  as.integer(1),  		# NMISSING
                  as.double(rep(missval,NQ)),	# KMISS
                  as.integer(minscale),		# MINSCALE
                  as.integer(rep(1,N)),		# MID
                  as.double(rawdata),		# KISSUE
                  fits = double(7*dims), 		# FITS
                  psimatrix = double(N*((dims*(dims+1))/2)+2*N*dims),		# PSIMATRIX
                  wmatrix = double((NQ)*((dims*(dims+1))/2)+2*(NQ)*dims),		# WMATRIX
                  lresp = integer(N+NQ),						# LRESPONDENTS
                  lmark = integer(N),		# LMARK
                  fits2 = double(6),		# FITS2
                  exitstatus = integer(1))	# EXITSTATUS

  if (res$exitstatus != 1) stop("\n\n\t====== Blackbox-Transpose did not execute properly ======\n\n")

  stimuli <- vector("list", dims)
  start <- 1
  end <- 2*NQ
  for(i in 1:dims){
    stimuli[[i]] <- as.data.frame(matrix(round(res$wmatrix[start:end],digits=3),nrow=NQ,ncol=i+1,byrow=T))
    colnames(stimuli[[i]]) <- c(paste("coord",1:i,"D",sep=""),"R2")
    rownames(stimuli[[i]]) <- stimnames
    stimuli[[i]] <- cbind(N=res$lresp[(length(res$lresp)-NQ+1):length(res$lresp)],stimuli[[i]])
    start <- end + 1
    end <- start + (i+2)*NQ - 1
  }

  individuals <- vector("list", dims)
  start <- 1
  end <- 3*N
  for(i in 1:dims){
    individuals[[i]] <- as.data.frame(matrix(round(res$psimatrix[start:end],digits=3),nrow=N,ncol=i+2,byrow=T))
    colnames(individuals[[i]]) <- c("c",paste("w",1:i,sep=""),"R2")
    if(!is.null(rownames(data))) rownames(individuals[[i]]) <- rownames(data)
    start <- end + 1
    end <- start + (i+3)*N - 1
    individuals[[i]][!res$lmark,] <- NA
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

  class(result) <- c("blackbt")
  if (verbose) cat("\n\n\tBlackbox-Transpose estimation completed successfully.\n\n")
  result
}

#' Blackbox-Transpose Summary
#'
#' \code{summary.blackbt} reads an \code{blackbt} object and prints a summary.
#'
#'
#' @param object a \code{blackbt} output object.
#' @param ...  further arguments to \code{print}.
#' @return A summary of a \code{blackbt} object. For each dimension, reports
#' all stimuli with coordinates, individuals used for scaling, and fit. Also
#' summarizes number of rows, columns, total data entries, number of missing
#' entries, percent missing data, and sum of squares.
#' @seealso '\link{blackbox_transpose}', '\link{LC1980}',
#' '\link{plotcdf}', '\link{LC1980_bbt}'.
#' @keywords multivariate
#' @export
#' @method summary blackbt
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' LCdat=LC1980[,-1]	#Dump the column of self-placements
#'
#'
#' ### This command conducts estimates, which we instead load using data()
#'
#' #LC1980_bbt <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
#' data(LC1980_bbt)
#'
#' plot(LC1980_bbt)
#' par(ask=TRUE)
#' plotcdf(LC1980_bbt)
#' summary(LC1980_bbt)
#'
summary.blackbt <- function(object, ...){

  x <- object
  if(!inherits(x, "blackbt"))  stop("Input is not of class blackbt.")
  cat("\n\nSUMMARY OF BLACKBOX TRANSPOSE OBJECT")
  cat("\n----------------------------------")
  for(i in 1:x$dims){
    cat("\n")
    print(x$stimuli[[i]],...)
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

#' Predict method of blackbt objects
#'
#' \code{predict.blackbt} reads an \code{blackbt} object and uses the estimates
#' to generate a matrix of predicted values.
#'
#'
#' @param object A \code{blackbox} output object.
#' @param dims Number of dimensions used in prediction. Must be equal to or
#' less than number of dimensions used in estimation.
#' @param ...  Ignored.
#' @return A matrix of predicted values generated from the parameters estimated
#' from a \code{blackbt} object.
#' @seealso '\link{blackbox_transpose}', '\link{LC1980}', '\link{LC1980_bbt}'
#' @keywords multivariate
#' @method predict blackbt
#' @export
#' @examples
#'
#' ## Estimate blackbt object from example and call predict function
#' data(LC1980)
#' data(LC1980_bbt)
#' prediction <- predict(LC1980_bbt, dims=2)
#'
#' ## Examine predicted vs. observed values for first 10 respondents
#' ## First column of LC1980 are self-placements, which are excluded
#' LC1980[1:10,-1]
#' prediction[1:10,]
#'
#' ## Check correlation across all predicted vs. observed, excluding missing values
#' prediction[which(LC1980[,-1] %in% c(0,8,9))] <- NA
#' cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
#'
predict.blackbt <- function(object, dims=1, ...){

  ## Error catch for dims
  if(!inherits(object, "blackbt"))  stop("Input is not of class blackbt.")
  if(dims < 1)  stop("dims must be great than 1")
  if(dims > object$dims)  stop(paste("dims must be equal or less than the number of estimate dimensions, which is", object$dims))

  W.hat<- as.matrix(object$individuals[[dims]][,paste("w", 1:dims, sep="")])
  Psi.hat <- as.matrix(object$stimuli[[dims]][,paste("coord", 1:dims, "D", sep="")])
  Jn <- matrix(rep(1, nrow(Psi.hat)),ncol=1)
  c <- matrix(object$individuals[[dims]][,"c"],nrow=1)

  ## In blackbox, Keith already postmultiplied by sqrt(singular), so no need to do it here
  ## This is not the case in blackbox_transpose, which needs to be multiplied by sqrt(singular)
  singular <- sqrt(diag(object$fits$singular[1:dims]))
  X.hat <-  Psi.hat %*% singular %*% t(W.hat) + Jn %*% c
  X.hat <- t(X.hat)

  colnames(X.hat) <- rownames(object$stimuli[[dims]])
  rownames(X.hat) <- rownames(object$individuals[[dims]])

  return(X.hat)
}

#' Blackbox Transpose Coordinate Distribution Plot
#'
#' \code{plot.blackbt} reads an \code{blackbt} object and plots the probability
#' distribution of the respondents and stimuli.
#'
#'
#' @param x an \code{blackbt} output object.
#' @param xlim vector of length 2, fed to the \code{plot} function as the
#' \code{xlim} argument, which sets the minimum and maximum range of the
#' x-axis.
#' @param ...  other arguments to \code{plot}.
#' @return A plot of the probability distribution of the respondent ideal
#' points, along with the locations of the stimuli.
#' @seealso '\link{blackbox_transpose}', '\link{LC1980}',
#' '\link{plotcdf}',  '\link{LC1980_bbt}'.
#' @keywords multivariate
#' @export
#' @method plot blackbt
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' LCdat=LC1980[,-1]	#Dump the column of self-placements
#'
#'
#' ### This command conducts estimates, which we instead load using data()
#'
#' #LC1980_bbt <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
#' data(LC1980_bbt)
#'
#'
#'
#' plot(LC1980_bbt)
#' par(ask=TRUE)
#' plotcdf(LC1980_bbt)
#' summary(LC1980_bbt)
#'
#'
plot.blackbt <- function(x, xlim=c(-1,1), ...){

  if(!inherits(x, "blackbt"))  stop("Input is not of class blackbt.")

  colchoice <- rep(palette(),3)
  dens <- density(x$individuals[[1]]$w1,na.rm=TRUE)
  ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

  plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
       ylab="Density", main="Stimuli and Population Distribution", lab=c(5,5,7),
       col="blue", cex.main=1.2, cex.lab=1.2)

  for(i in 1:nrow(x$stimuli[[1]])){
    arrows(x$stimuli[[1]]$coord1D[i], 0.08, x$stimuli[[1]]$coord1D[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
    text(x=x$stimuli[[1]]$coord1D[i], y=0.1, rownames(x$stimuli[[1]])[i], srt=90, adj=0, font=2)
  }

  text(x=0.85*xlim[2],y=0.85*ymax,paste("N =",x$Ncol),cex=1.3)
}

#' @method plotcdf blackbt
#' @export
plotcdf.blackbt <- function(x, ..., align=NULL, xlim=c(-1.2,1)){

  if(!inherits(x, "blackbt"))  stop("Input is not of class blackbt.")

  colchoice <- rep(palette(),3)
  cdf <- ecdf(x$individuals[[1]]$w1)
  plot(cdf, lwd=2, xlim=xlim, lab=c(5,5,7), bty="n",
       cex.points=0.5, cex.main=1.2, cex.lab=1.2,
       xlab="Location", ylab="Cumulative Density",
       main="Stimuli and Population CDF", ...)

  abline(v=x$stimuli[[1]]$coord1D,lty=2,col='gray70')

  if(is.null(align)) align <- (min(cdf(x$stimuli[[1]]$coord1D))+xlim[1])/1.5
  for(i in 1:nrow(x$stimuli[[1]])){
    arrows(align, cdf(x$stimuli[[1]]$coord1D[i]), x$stimuli[[1]]$coord1D[i]-0.1, cdf(x$stimuli[[1]]$coord1D[i]), length=0.1, angle=20, col=colchoice[i], lwd=2)
    text(x=align, y=cdf(x$stimuli[[1]]$coord1D[i]), rownames(x$stimuli[[1]])[i], adj=1, font=2)
  }
}

#' Bootstrap of Blackbox Transpose Scaling
#'
#' \code{boot_blackbt} is a function automates the non-parametric bootstrapping
#' of \code{blackbox_transpose}.  The original function takes a matrix of
#' perceptual data, such as liberal-conservative rankings of various stimuli,
#' and recovers the true location of those stimuli in a spatial model. The
#' bootstrap simply applies this estimator across multiple resampled data sets
#' and stores the results of each iteration in a matrix.  These results can be
#' used to estimate uncertainty for various parameters of interest, and can be
#' plotted using the \code{plot.boot_blackbt} function.
#'
#'
#' @param data matrix of numeric values, containing the perceptual data.
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
#' @param dims integer, specifies the number of dimensions to be estimated.
#' @param dim.extract integer, specifies which dimension to extract results for
#' the bootstrap from.
#' @param minscale integer, specifies the minimum number of responses a
#' respondent needs needs to provide to be used in the scaling.
#' @param iter integer, number of iterations the bootstrap should run for.
#' @return An object of class \code{boot_blackbt}. This is simply a matrix of
#' dimensions iter x number of stimuli.  Each row stores the estimated stimuli
#' locations for each iteration.
#' @seealso \link{blackbox_transpose}.
#' @keywords multivariate
#' @export
#' @examples
#'
#'
#' data(LC1980)
#'
#' data=LC1980[,-1]
#'
#' # Not run to save time, but loaded object is the output
#'
#' # bootbbt <- boot_blackbt(data, missing=c(0,8,9), dims=1, minscale=8, iter=10)
#'
#' data("bootbbt")
#'
#'
#' plot(bootbbt)
#'
#'
#'
boot_blackbt <- function(data, missing=NULL, dims=1, dim.extract=dims, minscale, iter=100){

  if(dim.extract > dims) stop("dim.extract must be less than dims")

  original <- blackbox_transpose(data=data, missing=missing, dims=dims, minscale=minscale, verbose=FALSE)
  samples <- matrix(NA, nrow=iter, ncol=nrow(original$stimuli[[dims]]))
  getdim <- paste("coord", dim.extract, "D",sep="")
  samples[1,] <- unlist(original$stimuli[[dims]][getdim])
  colnames(samples) <- rownames(original$stimuli[[dims]])
  rownames(samples) <- 1:iter

  for(i in 2:iter){
    tmpdat <- data[sample(1:nrow(data),nrow(data), replace=TRUE),]
    rownames(tmpdat) <- NULL
    tmpres <- blackbox_transpose(tmpdat, missing=missing, dims=dims,
                                 minscale=minscale, verbose=FALSE)
    samples[i,] <- unlist(tmpres$stimuli[[dims]][getdim])

    if(i %% 10 ==0){
      cat("\n\t\tIteration", i, "complete...")
      flush.console()
    }
  }

  class(samples) <- "boot_blackbt"
  samples <- sign(samples[,1])*samples*-1	# For polarity consistency
  return(samples)
}

#' Bootstrapped Blackbox Transpose Stimulus Plots
#'
#' \code{plot.boot_blackbt} reads an \code{boot_blackbt} object and plots a
#' dotchart of the stimuli with estimated confidence intervals.
#'
#'
#' @param x an \code{boot_blackbt} output object.
#' @param ...  other arguments to \code{plot}.
#' @return
#'
#' A dotchart of estimated stimulus locations, with 95 percent confidence
#' intervals. Point estimates are estimates from the original data set.
#' @seealso \link{blackbox_transpose}, \link{boot_blackbt}.
#' @keywords multivariate
#' @export
#' @method plot boot_blackbt
#' @examples
#'
#'
#'
#' data(LC1980)
#' data=LC1980[,-1]
#'
#' # Not run to save time, but loaded object is the output
#' # bootbbt <- boot_blackbt(data, missing=c(0,8,9), dims=1, minscale=8, iter=10)
#' data("bootbbt")
#'
#' plot(bootbbt)
#'
#'
plot.boot_blackbt <- function(x, ...){

  if(!inherits(x, "boot_blackbt"))  stop("Input is not of class boot_blackbt.")

  x <- x[,order(colMeans(x))]
  positions <- colMeans(x)
  xrange <- ceiling(max(abs(positions))*10)/10
  plot(positions, 1:ncol(x), type="n", xlim=c(-xrange, xrange), yaxt="n",bty="n",
       ylim=c(0.8,ncol(x)),xlab="Blackbox Transpose Score",ylab="",cex.lab=1.2, ...)
  bars <- apply(x,2,quantile,c(0.025,0.975))
  for(i in 1:ncol(bars)) segments(bars[1,i], i, bars[2,i],col="grey",lwd=2)
  points(positions, 1:ncol(bars), pch=20, cex=0.6)
  mtext(colnames(x), side=2, at=1:ncol(x),las=1,adj=1,cex=0.8, font=3)

}


