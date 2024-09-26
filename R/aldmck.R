#' Aldrich-McKelvey Scaling
#'
#' \code{aldmck} is a function that takes a matrix of perceptual data, such as
#' liberal-conservative rankings of various stimuli, and recovers the true
#' location of those stimuli in a spatial model. It differs from procedures
#' such as \code{wnominate}, which instead use preference data to estimate
#' candidate and citizen positions. The procedure here, developed by John
#' Aldrich and Richard McKelvey in 1977, is restricted to estimating data with
#' no missing values and only in one dimension.  Please refer to the
#' \code{blackbox} and \code{blackbox_transpose} functions in this package for
#' procedures that accomodate missing data and multidimensionality estimates.
#'
#'
#' @param data matrix of numeric values, containing the perceptual data.
#' Respondents should be organized on rows, and stimuli on columns. It is
#' helpful, though not necessary, to include row names and column names.
#' @param respondent integer, specifies the column in the data matrix of the
#' stimuli that contains the respondent's self-placement on the scale. Setting
#' respondent = 0 specifies that the self-placement data is not available.
#' Self-placement data is not required to estimate the locations of the
#' stimuli, but is required if recovery of the respondent ideal points, or
#' distortion parameters is desired. Note that no distortion parameters are
#' estimated in AM without self-placements because they are not needed, see
#' equation (24) in Aldrich and McKelvey (1977) for proof.
#' @param missing vector or matrix of numeric values, sets the missing values
#' for the data.  NA values are always treated as missing regardless of what is
#' set here.  Observations with missing data are discarded before analysis.  If
#' input is a vector, then the vector is assumed to contain the missing value
#' codes for all the data.  If the input is a matrix, it must be of dimension p
#' x q, where p is the maximum number of missing values and q is the number of
#' columns in the data.  Each column of the inputted matrix then specifies the
#' missing data values for the respective variables in data. If null (default),
#' no missing values are in the data other than the standard NA value.
#' @param polarity integer, specifies the column in the data matrix of the
#' stimuli that is to be set on the left side (generally this means a liberal)
#' @param verbose logical, indicates whether \code{aldmck} should print out
#' detailed output when scaling the data.
#' @return An object of class \code{aldmck}.
#'
#' \item{legislators}{ vector, containing the recovered locations of the
#' stimuli. The names of the stimuli are attached if provided as column names
#' in the argument \code{data}, otherwise they are generated sequentiall as
#' 'stim1', 'stim2', etc. }
#'
#' \item{respondents}{ matrix, containing the information estimated for each
#' respondent.  Observations which were discarded in the estimation for missing
#' data purposes have been NA'd out: \itemize{ \item\code{intercept}Intercept
#' of perceptual distortion for respondent.  \item\code{weight}Weight of
#' perceptual distortion for respondent.  \item\code{idealpt}Estimated location
#' of the respondent. Note that these positions are still calculated for
#' individuals with negative weights, so these may need to be discarded. Note
#' that this will not be calculated if self-placements are not provided in the
#' data.  \item\code{selfplace}The self-reported location of the individual,
#' copied from the \code{data} argument if \code{respondent} is not set to 0.
#' \item\code{polinfo}Estimated political information of respondent, calculated
#' as the correlation between the true and reported stimulus locations.  The
#' validation of this measure is provided in the article by Palfrey and Poole
#' in the references.  Note that this measure is included even for respondents
#' that were not used in the estimation. Individuals with negative weights have
#' also been assigned a political information score of 0, rather than negative
#' scores.  } } \item{eigenvalues}{ A vector of the eigenvalues from the
#' estimation.} \item{AMfit}{ Ratio of overall variance to perceptions in
#' scaled data divided by average variance. This measure of fit, described by
#' Aldrich and McKelvey, measures the amount of reduction of the variance of
#' the scaled over unscaled data.} \item{N}{ Number of respondents used in the
#' estimation (i.e. had no missing data)} \item{N.neg}{ Number of cases with
#' negative weights. Only calculated if respondent self-placements are
#' specified, will equal 0 if not.} \item{N.pos}{ Number of cases with positive
#' weights. Only calculated if respondent self-placements are specified, will
#' equal 0 if not.}
#' @seealso '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.aldmck}',
#' '\link{plotcdf}'.
#' @references
#'
#' Keith Poole, Jeffrey Lewis, Howard Rosenthal, James Lo, Royce Carroll (2016)
#' ``Recovering a Basic Space from Issue Scales in R.'' Journal of Statistical
#' Software. 69(7), 1--21. doi:10.18637/jss.v069.i07
#'
#' John H. Aldrich and Richard D. McKelvey (1977) ``A Method of Scaling with
#' Applications to the 1968 and 1972 Presidential Elections.'' American
#' Political Science Review. 71(1), 111-130.
#'
#' Thomas R. Palfrey and Keith T. Poole (1987) ``The Relationship between
#' Information, Ideology, and Voting Behavior.'' American Journal of Political
#' Science. 31(3), 511-530.
#'
#' Keith Poole. \url{https://voteview.com}
#' @export
#' @keywords multivariate
#' @examples
#'
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' plot(result)
#'
#'
aldmck <- function(data, respondent = 0, missing=NULL, polarity, verbose=FALSE) {

  ### Error check each argument ###
  if(inherits(data, "data.frame"))  data <- as.matrix(data)
  if(!inherits(data, "matrix"))  stop("Data is not a matrix, or data frame.")
  if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric() or mode().")

  if(mode(respondent) != "numeric")  stop("Respondent is not specified as an integer.")
  if(respondent > ncol(data))  stop("Respondent set to a column greater than number of columns in data.")
  if(respondent < 0)  stop("Respondent cannot be negative, set respondent=0 is self identification is unavailable.")

  if(mode(polarity) != "numeric")  stop("Polarity is not specified as an integer.")
  if(polarity == respondent) stop("Self-placements cannot also be set as polarity.")
  if(polarity > ncol(data))  stop("Polarity set to a column greater than number of columns in data.")
  if(polarity < 0)  stop("Polarity cannot be negative.")

  if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
  if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")

  if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")

  ### Format the data input ###
  N <- nrow(data)
  NQ <- ncol(data) - 1
  NRESP <- respondent
  if(respondent==0){
    NQ <- ncol(data)
    NRESP  <- ncol(data) + 1
  }

  if(is.vector(missing))	data[data %in% missing] <- NA
  if(is.matrix(missing))	for(i in 1:ncol(data))	data[data[,i] %in% missing[,i],i] <- NA
  missval <- max(data,na.rm=TRUE) + 1	#code for missing data
  if(respondent==0) data <- cbind(data, fakedata = rep(missval,nrow(data)))
  rawdata <- as.integer(t(data))
  rawdata[is.na(rawdata)] <- missval

  ##Longer output
  deleted <- sum(is.na(apply(data,1,sum)))
  stimnames <- colnames(data)
  if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
  if(verbose){

    cat("\n\n\tBeginning Aldrich-McKelvey Scaling...")
    if(respondent!=0) cat("\n\n\t\tColumn '",stimnames[respondent],"' is set as the self placement.", sep="")
    if(respondent==0) cat("\n\n\t\tSelf-placements have not been provided by the user.")

    cat("\n\t\tColumn '",stimnames[polarity],"' is set as the left-leaning stimulus.\n\t\t", sep="")
    cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
    cat(NQ, "stimuli have been provided.")
  }
  if ((respondent != 0) & (respondent < polarity))   polarity = polarity - 1


  ### Send the data to Fortran ###
  res <- .Fortran("mckalnew",
                  as.integer(N), 			# NRESPONDENTS
                  as.integer(NQ+1),	 		# NISSUES
                  as.integer(NRESP),  		# NSELFPOS
                  as.integer(1),			# NMISSING (maximum number of missing valaues)
                  as.numeric(rep(missval,NQ+1)),	# KMISS, input number of missing values in each stimuli, length(NMISS*(NQ+1))), need to be double later
                  as.integer(polarity),       	# POLARITY (input, which is on left)
                  as.integer(rep(1,N)),		# MID, input, ID number of individual,length(nrow(data))
                  as.numeric(rawdata),		# KISSUE, data matrix, needs to be double later, length(nrow(data)*(NQ+1)))
                  #              as.character("a"),		# CAND, input, names of stimuli each column,length(NQ+1)
                  fits = double(5),			# FITS (output), AMfit, R^2, numberscaled, negweights, posweights
                  psimatrix = double(N*4),		# PSIMATRIX, intercept C, W weight, scaleposition, R^2
                  stimcoord = double(NQ),		# STIMCOORDS, stimulus coordinates
                  eigenvalues = double(NQ),		# EIGENVALUES
                  exitstatus = integer(1))		# EXITSTATUS
  if(res$exitstatus!=1) stop("\n\n\t====== Aldrich-McKelvey did not execute properly ======\n\n")

  ### Format the output ###
  stimnames <- colnames(data)
  if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
  stimuli <- as.numeric(res$stimcoord)
  names(stimuli) <- stimnames[-NRESP]

  respondents <- matrix(res$psimatrix, ncol=4, byrow=T)
  respondents[respondents==0] <- NA
  midnames <- as.character(unlist(rownames(data)))
  if(is.null(midnames)) midnames <- paste("resp", 1:N, sep="")
  rownames(respondents) <- midnames
  colnames(respondents) <- c("intercept", "weight", "idealpt", "R2")

  ## For JSS edit, R2 removed
  respondents <- respondents[,1:3]

  if (respondent == 0) {
    selfplace <- rep(NA, N)
    data <- data[,-ncol(data)]
  }    else {
    selfplace <- data[, respondent]
    data <- data[, -respondent]
  }

  polinfo <- suppressWarnings(apply(data,1,cor,y=stimuli,use="pairwise.complete.obs"))
  polinfo[which(respondents[,"weight"] < 0)] <- 0
  respondents <- cbind(respondents, selfplace, polinfo)


  result <- list(stimuli = stimuli, respondents = as.data.frame(respondents),
                 eigenvalues = as.numeric(res$eigenvalues),
                 AMfit = as.numeric(res$fits[1]),
                 #			R2 = as.numeric(res$fits[2]),
                 #			N = nrow(data)-deleted, N.neg = as.integer(res$fits[5]),
                 N = as.integer(res$fits[4])+as.integer(res$fits[5]), N.neg = as.integer(res$fits[5]),
                 N.pos = as.integer(res$fits[4]))

  class(result) <- c("aldmck")
  if(verbose) cat("\n\n\tAldrich-McKelvey estimation completed successfully.\n\n")
  result
}

#' Aldrich-McKelvey Summary
#'
#' \code{summary.aldmck} reads an \code{aldmck} object and prints a summary.
#'
#'
#' @param object an \code{aldmck} output object.
#' @param ...  further arguments to \code{print}.
#' @return A summary of an \code{aldmck} object. Reports number of stimuli,
#' respondents scaled, number of respondents with positive and negative
#' weights, R-squared, Reudction of normalized variance of perceptions, and
#' stimuli locations.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{plot.aldmck}',
#' '\link{plotcdf}'.
#' @keywords multivariate
#' @export
#' @method summary aldmck
#' @examples
#'
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' par(ask=TRUE)
#' plot.AM(result,xlim=c(-1.5,1.5))
#' plotcdf(result)
#'
summary.aldmck <- function(object, ...){

  x<-object
  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")


  cat("\n\nSUMMARY OF ALDRICH-MCKELVEY OBJECT")
  cat("\n----------------------------------")
  cat("\n\nNumber of Stimuli:", length(x$stimuli))
  cat("\nNumber of Respondents Scaled:", x$N)
  if(sum(!is.na(x$respondents[,"selfplace"])) != 0){
    cat("\nNumber of Respondents (Positive Weights):", x$N.pos)
    cat("\nNumber of Respondents (Negative Weights):", x$N.neg)
    #	cat("\n\nR-Squared:", round(x$R2, digits=2)
  }
  cat("\nReduction of normalized variance of perceptions:", round(x$AMfit, digits=2), "\n\n")

  final <- matrix(round(sort(x$stimuli),3),ncol=1)
  rownames(final) <- names(sort(x$stimuli))
  colnames(final) <- c("Location")
  print(final,...)
  cat("\n\n")
}

#' Predict method of aldmck objects
#'
#' \code{predict.aldmck} reads an \code{aldmck} object and uses the estimates
#' to generate a matrix of predicted values.
#'
#'
#' @param object A \code{aldmck} output object.
#' @param caliper Caliper tolerance. Any individuals with estimated weights
#' lower than this value are NA'd out for prediction. Since predictions are
#' made by dividing observed values by estimating weights, very small weights
#' will grossly inflate the magnitude of predicted values and lead to extreme
#' predictions.
#' @param ...  Ignored.
#' @return A matrix of predicted values generated from the parameters estimated
#' from a \code{aldmck} object.
#' @seealso '\link{aldmck}', '\link{LC1980}'
#' @keywords multivariate
#' @export
#' @method predict aldmck
#' @examples
#'
#' ## Estimate an aldmck object from example and call predict function
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' prediction <- predict(result)
#'
#' ## Examine predicted vs. observed values for first 10 respondents
#' ## Note some observations are NA'd in prediction matrix from caliper
#' ## First column of LC1980 are self-placements, which are excluded
#' LC1980[1:10,-1]
#' prediction[1:10,]
#'
#' ## Check correlation across all predicted vs. observed, excluding missing values
#' prediction[which(LC1980[,-1] %in% c(0,8,9))] <- NA
#' cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
#'
predict.aldmck <- function(object, caliper=0.2, ...){

  if(!inherits(object, "aldmck"))  stop("Input is not of class aldmck.")
  if(caliper < 0) stop("Caliper must be positive.")

  N <-nrow(object$respondents)	#number of respondents
  j <- length(object$stimuli)	#number of stimuli
  Y <- matrix(rep(object$stimuli,N), nrow=N, ncol=j, byrow=TRUE)
  c <- matrix(object$respondent[,"intercept"],ncol=1) %*% matrix(rep(1,j),nrow=1)
  w <- matrix(object$respondent[,"weight"],ncol=1) %*% matrix(rep(1,j),nrow=1)
  X.hat <- (Y - c)/w
  dumpthese <- which(abs(object$respondent[,"weight"]) < caliper)
  X.hat[dumpthese,] <- NA

  colnames(X.hat) <- names(object$stimuli)
  rownames(X.hat) <- rownames(object$respondents)

  return(X.hat)
}


#' Aldrich-McKelvey Coordinate Distribution Plot
#'
#' \code{plot.aldmck} reads an \code{aldmck} object and plots the probability
#' distribution of the respondents and stimuli.
#'
#' @param x an \code{aldmck} output object.
#' @param ...  Other arguments to \code{plot}.
#' @return A plot of the probability distribution of the respondent ideal
#' points, along with the locations of the stimuli. If no self-placements were
#' specified during estimation, no graphical plots will appear.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}',
#' '\link{plot.AM}', '\link{plotcdf}'
#' '\link{plot.aldmck_negative}','\link{plot.aldmck_positive}'.
#' @keywords multivariate
#' @export
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' plot(result)
#'
plot.aldmck <- function(x, ...){

  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")

  op <- par(no.readonly=TRUE)
  par(mfrow=c(2,2))
  plot.AM(x,...)
  plotcdf(x,...)
  plot.aldmck_positive(x,...)
  plot.aldmck_negative(x,...)
  suppressWarnings(par(op))

}



#' Aldrich-McKelvey Coordinate Distribution Plot
#'
#' \code{plot.AM} reads an \code{aldmck} object and plots the probability
#' distribution of the respondents and stimuli.
#'
#'
#' @param x an \code{aldmck} output object.
#' @param xlim vector of length 2, fed to the \code{plot} function as the
#' \code{xlim} argument, which sets the minimum and maximum range of the
#' x-axis.
#' @param ...  other arguments to \code{plot}.
#' @return A plot of the probability distribution of the respondent ideal
#' points, along with the locations of the stimuli. If no self-placements were
#' specified during estimation, no graphical plots will appear.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}',
#' '\link{plotcdf}', '\link{plot.aldmck}'
#' @keywords multivariate
#' @rawNamespace export(plot.AM)
#' @exportS3Method NULL
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' par(ask=TRUE)
#' plot.AM(result,xlim=c(-1.5,1.5))
#' plotcdf(result)
#'
plot.AM <- function(x, xlim=c(-2,2), ...){

  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")

  if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
    plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
    text(0.5,0.5,"No\nself\nplacement",cex=3)
    return()
  }

  colchoice <- rep(palette(),3)
  dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
  ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

  plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
       ylab="Density", main="Stimuli and Population Distribution", lab=c(5,5,7),
       col="blue", cex.main=1.2, cex.lab=1.2, ...)

  for(i in 1:length(x$stimuli)){
    arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
    text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
  }

  text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N),cex=1.3)
}




#' Aldrich-McKelvey Positive Coordinate Distribution Plot
#'
#' \code{plot.aldmck_positive} reads an \code{aldmck} object and plots the
#' probability distribution of the respondents and stimuli with positive
#' weights.
#'
#'
#' @param x an \code{aldmck} output object.
#' @param xlim vector of length 2, fed to the \code{plot} function as the
#' \code{xlim} argument, which sets the minimum and maximum range of the
#' x-axis.
#' @param ...  other arguments to \code{plot}.
#' @return A plot of the probability distribution of the respondent ideal
#' points, along with the locations of the stimuli. If no weights exist because
#' respondent self-placements are not specified, a plot indicating this in text
#' is given.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}',
#' '\link{plotcdf}', '\link{plot.aldmck}'
#' @keywords multivariate
#' @rawNamespace export(plot.aldmck_positive)
#' @exportS3Method NULL
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' plot.aldmck_positive(result,xlim=c(-1.5,1.5))
#'
plot.aldmck_positive <- function(x, xlim=c(-2,2), ...){

  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")

  if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
    plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
    text(0.5,0.5,"No\nweights",cex=3)
    return()
  }

  x$respondents <- x$respondents[x$respondents[,"weight"]>0,]

  colchoice <- rep(palette(),3)
  dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
  ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

  plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
       ylab="Density", main="Positive Weights", lab=c(5,5,7),
       col="blue", cex.main=1.2, cex.lab=1.2, ...)

  for(i in 1:length(x$stimuli)){
    arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
    text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
  }

  text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.pos),cex=1.3)
}



#' Aldrich-McKelvey Negative Coordinate Distribution Plot
#'
#' \code{plot.aldmck_negative} reads an \code{aldmck} object and plots the
#' probability distribution of the respondents and stimuli with negative
#' weights.
#'
#'
#' @param x an \code{aldmck} output object.
#' @param xlim vector of length 2, fed to the \code{plot} function as the
#' \code{xlim} argument, which sets the minimum and maximum range of the
#' x-axis.
#' @param ...  other arguments to \code{plot}.
#' @return A plot of the probability distribution of the respondent ideal
#' points, along with the locations of the stimuli. If no negative weights
#' exist, either because respondent self-placements are not specified, or
#' because all weights are positive, a plot indicating this in text is given.
#' @seealso '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}',
#' '\link{plotcdf}', '\link{plot.aldmck}'
#' @keywords multivariate
#' @rawNamespace export(plot.aldmck_negative)
#' @exportS3Method NULL
#' @examples
#'
#' ### Loads and scales the Liberal-Conservative scales from the 1980 NES.
#' data(LC1980)
#' result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
#' summary(result)
#' plot.aldmck_negative(result,xlim=c(-1.5,1.5))
#'
plot.aldmck_negative <- function(x, xlim=c(-2,2), ...){

  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")

  if(sum(x$respondents[,"weight"]<0,na.rm=TRUE)==0){
    plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
    text(0.5,0.5,"No\nnegative\nweights",cex=2.5)
    return()
  }

  x$respondents <- x$respondents[x$respondents[,"weight"]<0,]
  dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
  ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

  if(x$N.neg>50){
    colchoice <- rep(palette(),3)
    dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
    ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

    plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
         ylab="Density", main="Negative Weights", lab=c(5,5,7),
         col="blue", cex.main=1.2, cex.lab=1.2, ...)

    for(i in 1:length(x$stimuli)){
      arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
      text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
    }
    text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.neg),cex=1.3)
  } #end if N>50

  if(x$N.neg<51){
    hist(x$respondents[,"idealpt"],xlim=c(-2,2),ylim=c(0,ymax),breaks=30,col="lightgrey",freq=F,
         xlab="Location",ylab="Density", main="Stimuli and Population Distribution: Negative Weights", cex.main=1.2, cex.lab=1.2,...)
    text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.neg),cex=1.3)

  } #end if N<51
}



#' @method plotcdf aldmck
#' @export
plotcdf.aldmck <- function(x, ..., align=NULL, xlim=c(-2,2)){

  if(!inherits(x, "aldmck"))  stop("Input is not of class aldmck.")

  if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
    plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
    text(0.5,0.5,"No\nself\nplacement",cex=3)
    return()
  }

  colchoice <- rep(palette(),3)
  cdf <- ecdf(x$respondents[,"idealpt"])
  plot(cdf, lwd=2, xlim=xlim, lab=c(5,5,7), bty="n",
       cex.points=0.5, cex.main=1.2, cex.lab=1.2,
       xlab="Location", ylab="Cumulative Density",
       main="Stimuli and Population CDF", ...)

  abline(v=x$stimuli,lty=2,col='gray70')

  if(is.null(align)) align <- cdf(x$stimuli[1]) - 1.5
  for(i in 1:length(x$stimuli)){
    arrows(align + 0.1, cdf(x$stimuli[i]), x$stimuli[i]-0.1, cdf(x$stimuli[i]), length=0.1, angle=20, col=colchoice[i], lwd=2)
    text(x=align, y=cdf(x$stimuli[i]), names(x$stimuli)[i], adj=1, font=2)
  }
}


#' Bootstrap of Aldrich-McKelvey Scaling
#'
#' \code{boot_aldmck} is a function automates the non-parametric bootstrapping
#' of \code{aldmck}.  The original function takes a matrix of perceptual data,
#' such as liberal-conservative rankings of various stimuli, and recovers the
#' true location of those stimuli in a spatial model. The bootstrap simply
#' applies this estimator across multiple resampled data sets and stores the
#' results of each iteration in a matrix.  These results can be used to
#' estimate uncertainty for various parameters of interest, and can be plotted
#' using the \code{plot.boot_aldmck} function.
#'
#'
#' @param data matrix of numeric values, containing the perceptual data.
#' Respondents should be organized on rows, and stimuli on columns. It is
#' helpful, though not necessary, to include row names and column names.
#' @param respondent integer, specifies the column in the data matrix of the
#' stimuli that contains the respondent's self-placement on the scale. Setting
#' respondent = 0 specifies that the self-placement data is not available.
#' Self-placement data is not required to estimate the locations of the
#' stimuli, but is required if recovery of the respondent ideal points, or
#' distortion parameters is desired. Note that no distortion parameters are
#' estimated in AM without self-placements because they are not needed, see
#' equation (24) in Aldrich and McKelvey (1977) for proof.
#' @param missing vector or matrix of numeric values, sets the missing values
#' for the data.  NA values are always treated as missing regardless of what is
#' set here.  Observations with missing data are discarded before analysis.  If
#' input is a vector, then the vector is assumed to contain the missing value
#' codes for all the data.  If the input is a matrix, it must be of dimension p
#' x q, where p is the maximum number of missing values and q is the number of
#' columns in the data.  Each column of the inputted matrix then specifies the
#' missing data values for the respective variables in data. If null (default),
#' no missing values are in the data other than the standard NA value.
#' @param polarity integer, specifies the column in the data matrix of the
#' stimuli that is to be set on the left side (generally this means a liberal)
#' @param iter integer, is the number of iterations the bootstrap should run
#' for.
#' @return An object of class \code{boot_aldmck}. This is simply a matrix of
#' dimensions iter x number of stimuli.  Each row stores the estimated stimuli
#' locations for each iteration.
#' @seealso '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.aldmck}',
#' '\link{plotcdf}'.
#' @references
#'
#' Keith Poole, Jeffrey Lewis, Howard Rosenthal, James Lo, Royce Carroll (2016)
#' ``Recovering a Basic Space from Issue Scales in R.'' Journal of Statistical
#' Software. 69(7), 1--21. doi:10.18637/jss.v069.i07
#'
#' John H. Aldrich and Richard D. McKelvey (1977) ``A Method of Scaling with
#' Applications to the 1968 and 1972 Presidential Elections.'' American
#' Political Science Review. 71(1), 111-130.
#'
#' Thomas R. Palfrey and Keith T. Poole (1987) ``The Relationship between
#' Information, Ideology, and Voting Behavior.'' American Journal of Political
#' Science. 31(3), 511-530.
#'
#' Keith Poole. \url{https://voteview.com}
#' @keywords multivariate
#' @export
#' @examples
#'
#'
#' data(LC1980)
#'
#' result <- boot_aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), iter=30)
#' plot(result)
boot_aldmck <- function(data, respondent=0, missing=NULL, polarity, iter=100){

  original <- aldmck(data=data, respondent=respondent, polarity=polarity,
                     missing=missing,verbose=FALSE)
  samples <- matrix(NA, nrow=iter, ncol=length(original$stimuli))
  samples[1,] <- original$stimuli
  colnames(samples) <- names(original$stimuli)
  rownames(samples) <- 1:iter

  for(i in 2:iter){
    tmpdat <- data[sample(1:nrow(data),nrow(data), replace=TRUE),]
    samples[i,] <- aldmck(data=tmpdat, respondent=respondent,
                          polarity=polarity, missing=missing,verbose=FALSE)$stimuli
  }

  class(samples) <- "boot_aldmck"
  samples <- sign(samples[,1])*samples*-1	# For polarity consistency
  return(samples)
}



#' Bootstrapped Aldrich-McKelvey Stimulus Plots
#'
#' \code{plot.boot_aldmck} reads an \code{boot_aldmck} object and plots a
#' dotchart of the stimuli with estimated confidence intervals.
#'
#'
#' @param x an \code{boot_aldmck} output object.
#' @param ...  other arguments to \code{plot}.
#' @return
#'
#' A dotchart of estimated stimulus locations, with 95 percent confidence
#' intervals. Point estimates are estimates from the original data set.
#' @seealso \link{aldmck}, \link{boot_aldmck}.
#' @keywords multivariate
#' @export
#' @method plot boot_aldmck
#' @examples
#'
#'
#'
#' data(LC1980)
#' result <- boot_aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), iter=30)
#' plot(result)
#'
#'
plot.boot_aldmck <- function(x, ...){

  if(!inherits(x, "boot_aldmck"))  stop("Input is not of class boot_aldmck.")

  x <- x[,order(colMeans(x))]
  positions <- colMeans(x)
  xrange <- ceiling(max(abs(positions))*10)/10
  plot(positions, 1:ncol(x), type="n", xlim=c(-xrange, xrange), yaxt="n",bty="n",
       ylim=c(0.8,ncol(x)),xlab="AM position score",ylab="",cex.lab=1.2, ...)
  bars <- apply(x,2,quantile,c(0.025,0.975))
  for(i in 1:ncol(bars)) segments(bars[1,i], i, bars[2,i],col="grey",lwd=2)
  points(positions, 1:ncol(bars), pch=20, cex=0.6)
  mtext(colnames(x), side=2, at=1:ncol(x),las=1,adj=1,cex=0.8, font=3)

}
