#' Blackbox Transpose Bootstrap of 1980 Liberal-Conservative Scales.
#'
#' Output from 10 bootstrap trials of LC1980 data. Included to allow the
#' example to run sufficiently quickly to pass CRAN guidelines.
#' @format ## `bootbbt`
#' A data frame with 10 rows and 6 columns giving the bootstrap placements for Carter, Reagan, Kennedy, Anderson, Republicans and Democrats
#' @keywords datasets
#'
"bootbbt"

#' 2004 PELA Liberal-Conservative Scales.
#'
#' Liberal-Conservative 10-point scales from the University of Salamanca's
#' Parliamentary Elites of Latin America (PELA) survey. Stored as a matrix of
#' integers. The number 99 is a missing value. These data come from Sebastian
#' Saiegh and are used in the paper and book cited below.
#'
#' @format ## `colombia`
#' An integer matrix with 95 rows and 17 columns:
#' \describe{
#'    \item{id}{Respondent ID.}
#'    \item{party}{Respondent party.}
#'    \item{departam}{Respondent district.}
#'    \item{entrev}{Interviewer ID.}
#'    \item{pl_uribista}{Placement of Partido Liberal Uribista on 10 point scale.}
#'    \item{pl_oficial}{Placement of Partido Liberal Oficial on 10 point scale.}
#'    \item{conservator}{Placement of Partido Conservador on 10 point scale.}
#'    \item{polo}{Placement of Polo on 10 point scale.}
#'    \item{union_cristiana}{Placement of Union Cristiana on 10 point scale.}
#'    \item{salvation}{Placement of Salvacion on 10 point scale.}
#'    \item{urine}{Placement of Mr. Uribe on 10 point scale.}
#'    \item{antanas}{Placement of Mr. Antanas on 10 point scale.}
#'    \item{gomez}{Placement of Mr. Gomez on 10 point scale.}
#'    \item{garzon}{Placement of Garzon on 10 point scale.}
#'    \item{holgin}{Placement of Holguin on 10 point scale.}
#'    \item{rivera}{Placement of Rivera on 10 point scale.}
#'    \item{self}{Respondent self placement on 10 point scale.}
#' }
#' @source
#'
#' Sebastian Saiegh. 2009. Recovering a Basic Space from Elite Surveys:
#' Evidence from Latin America. Legislative Studies Quarterly. 34(1): 117-145. \doi{10.3162/036298009787500349}
#'
#' Sebastian Saiegh. 2011. Ruling By Statute: How Uncertainty and Vote-Buying
#' Shape Lawmaking. New York: Cambridge University Press. \doi{10.1017/CBO9780511842276}
#' @keywords datasets
"colombia"

#' Blackbox Estimate, 1980 NES Issue Scales.
#'
#' Blackbox estimates from issues scales from the 1980 National Election Study.
#'
#' @format ## `Issues1980_bb`
#' An object of class `blackbox` with the following elements:
#' \describe{
#'    \item{stimuli}{A list of length 3 (number of dimensions in the scaling analysis) where each element contains a data frame with the following variables:}
#'    \itemize{
#'      \item {N}Number of respondents who provided a response to this stimulus.
#'      \item{c}Stimulus intercept.
#'      \item{w1}Estimate of the stimulus weight on the first dimension. If viewing the results for a higher dimension, higher dimension results will appear as w2, w3, etc.
#'      \item{R2}The percent variance explained for the stimulus. This increases as more dimensions are estimated.
#'    }
#'    \item{individuals}{A list of length 3 (number of dimensions in the scaling analysis) where each element is a data frame presenting results from the analysis fot he corresponding dimension.  The datasets contain the following variables:}
#'    \itemize{
#'      \item{c1}Estimate of the individual intercept on the first dimension. If viewing the results for a higher dimension, higher dimension results will appear as c2, c3, etc.
#'    }
#'    \item{fits}{ A data frame of fit results, with elements listed as follows:}
#'    \itemize{
#'      \item{SSE}Sum of squared errors.
#'      \item{SSE.explained}Explained sum of squared error.
#'      \item{percent}Percentage of total variance explained.
#'      \item{SE}Standard error of the estimate, with formula provided on pg. 973 of the article cited below.
#'      \item{singular}Singluar value for the dimension.
#'    }
#'    \item{Nrow}{Numberof rows/stimuli.}
#'    \item{Ncol}{Number of columns used in estimation. This may differ from the data set due to columns discarded due to the `minscale` constraint.}
#'    \item{Ndata}{Total number of data entries.}
#'    \item{Nmiss}{Number of missing entries.}
#'    \item{SS_mean}{ Sum of squares grand mean.}
#'    \item{dims}{ Number of dimensions estimated.}
#' }
#' @source
#' American National Election Study. \url{https://electionstudies.org}
#' @keywords datasets
#' @docType data
#' @name Issues1980_bb
#' @usage data(Issues1980_bb)
"Issues1980_bb"




#' 1980 Issues Scales
#'
#' Issue scales from the 1980 National Election Study. The numbers 0, 8, and 9
#' are considered to be missing values, except for the two abortion scales,
#' where '7' is also a missing value.  Hence, it must be recoded as in the
#' example shown below before scaling. The data is used as an example for
#' blackbox().
#'
#'
#' @format ## `Issues1980`
#' A matrix, containing reported self-placements along various stimuli on a  7 point
#' Liberal-Conservative scales (with the exception of abortion scales, which are 4 point):
#' \describe{
#'    \item{libcon1}{Liberal-conservative self-placement on 7 point scale.}
#'    \item{defense}{Defense spending self-placement on 7 point scale.}
#'    \item{govserv}{Government service on 7 point scale.}
#'    \item{inflation}{Importance of inflation self-placement on 7 point scale.}
#'    \item{abortion1}{Attitude on abortion 4 point scale.}
#'    \item{taxcut}{Support for tax cut on 7 point scale.}
#'    \item{libcon2}{Liberal-conservative self-placement on 7 point scale.}
#'    \item{govhelpmin}{Government aid on 7 point scale.}
#'    \item{russia}{Attitude towards Russia on 7 point scale.}
#'    \item{womenrole}{Role of women on 7 point scale.}
#'    \item{govjobs}{Placement of Democrats on 7 point scale.}
#'    \item{equalrights}{Support for equal rights on 7 point scale.}
#'    \item{busing}{Opinion on busing on 7 point scale.}
#'    \item{abortion2}{Another attitude on abortion on 4 point scale.}
#' }
#' @source
#' American National Election Study. \url{https://electionstudies.org}
#' Also available from Keith Poole's website. \url{https://voteview.com}
#' @keywords datasets
"Issues1980"





#' Blackbox Transpose Estimate, 1980 Liberal-Conservative Scales.
#'
#' Blackbox-Transpose estimates from Liberal-Conservative 7-point scales from
#' the 1980 National Election Study. Estimates in 3 dimensions.
#'
#'
#' @format ## `LC1980_bbt`
#' An object of class `blackbt` with the following elements
#' \describe{
#'    \item{stimuli}{A list of length 3 (number of dimensions in the scaling analysis) where each element contains a data frame with the following variables:}
#'    \itemize{
#'      \item{N}Number of respondents who provided a response to this stimulus.
#'      \item{coord1D}Location of the stimulus in the first dimension. If viewing the results for a higher dimension, higher dimension results will appear as coord2D, coord3D, etc.
#'      \item{R2}The percent variance explained for the stimulus. This increases as more dimensions are estimated.
#'    }
#'    \item{individuals}{A list of length 3 (number of dimensions in the scaling analysis) where each element is a data frame presenting results from the analysis fot he corresponding dimension.  The datasets contain the following variables:}
#'    \itemize{
#'      \item{c}Estimate of the individual intercept.
#'      \item{w1}Estimate of the individual slope. If viewing the results for a higher dimension, higher dimension results will appear as w2, w3, etc.
#'      \item{R2}The percent variance explained for the respondent. This increases as more dimensions are estimated.
#'    }
#'    \item{fits}{ A data frame of fit results, with elements listed as follows:}
#'    \itemize{
#'      \item{SSE}Sum of squared errors.
#'      \item{SSE.explained}Explained sum of squared error.
#'      \item{percent}Percentage of total variance explained.
#'      \item{SE}Standard error of the estimate, with formula provided on pg. 973 of the article cited below.
#'      \item{singular}Singluar value for the dimension.
#'    }
#'    \item{Nrow}{Numberof rows/stimuli.}
#'    \item{Ncol}{Number of columns used in estimation. This may differ from the data set due to columns discarded due to the `minscale` constraint.}
#'    \item{Ndata}{Total number of data entries.}
#'    \item{Nmiss}{Number of missing entries.}
#'    \item{SS_mean}{ Sum of squares grand mean.}
#'    \item{dims}{ Number of dimensions estimated.}
#' }
#' @source
#' American National Election Study. \url{https://electionstudies.org}
#' @keywords datasets
"LC1980_bbt"




#' 1980 Liberal-Conservative Scales.
#'
#' Liberal-Conservative 7-point scales from the 1980 National Election Study.
#' Includes (in order) self-placement, and rankings of Carter, Reagan, Kennedy,
#' Anderson, Republican party, Democratic Party. Stored as a matrix of
#' integers. The numbers 0, 8, and 9 are considered to be missing values.
#'
#' @format ## `LC1980`
#' matrix, containing reported placements of various stimuli on a 7 point Liberal-Conservative scale:
#' \describe{
#' \item{Self}{Self-placement on 7 point scale.}
#' \item{Carter}{Placement of Carter on 7 point scale.}
#' \item{Reagan}{Placement of Reagan on 7 point scale.}
#' \item{Kennedy}{Placement of Kennedy on 7 point scale.}
#' \item{Anderson}{Placement of Anderson on 7 point scale.}
#' \item{Republicans}{Placement of Republicans on 7 point scale.}
#' \item{Democrats}{Placement of Democrats on 7 point scale.}
#' }
#' @source
#' American National Election Study. \url{https://electionstudies.org}
#' Also available from Keith Poole's website. \url{https://voteview.com}
#' @keywords datasets
"LC1980"
