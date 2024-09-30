#include <Rcpp.h>
using namespace Rcpp;
// blackboxf_wrapper.cpp
#include "blackboxf.h" // Header file where your C++ function is declared

// [[Rcpp::export]]
SEXP blackboxf_wrapper(SEXP NRESPONDENTS_SEXP, SEXP NISSUES_SEXP, SEXP NDIMENSIONS_SEXP,
                       SEXP NMISSING_SEXP, SEXP KMISS_SEXP, SEXP MINSCALE_SEXP, SEXP MID_SEXP,
                       SEXP KISSUE_SEXP, SEXP FITS_SEXP, SEXP PSIMATRIX_SEXP, SEXP WMATRIX_SEXP,
                       SEXP LRESPONDENTS_SEXP, SEXP LMARK_SEXP, SEXP FITS2_SEXP, SEXP EXITSTATUS_SEXP)
{
  // Convert SEXP to C++ types
  int NRESPONDENTS = as<int>(NRESPONDENTS_SEXP);
  int NISSUES = as<int>(NISSUES_SEXP);
  int NDIMENSIONS = as<int>(NDIMENSIONS_SEXP);
  int NMISSING = as<int>(NMISSING_SEXP);
  NumericMatrix KMISS = as<NumericMatrix>(KMISS_SEXP);
  int MINSCALE = as<int>(MINSCALE_SEXP);
  IntegerVector MID = as<IntegerVector>(MID_SEXP);
  NumericMatrix KISSUE = as<NumericMatrix>(KISSUE_SEXP);
  NumericVector FITS = as<NumericVector>(FITS_SEXP);
  NumericMatrix PSIMATRIX = as<NumericMatrix>(PSIMATRIX_SEXP);
  NumericMatrix WMATRIX = as<NumericMatrix>(WMATRIX_SEXP);
  IntegerVector LRESPONDENTS = as<IntegerVector>(LRESPONDENTS_SEXP);
  IntegerVector LMARK = as<IntegerVector>(LMARK_SEXP);
  NumericVector FITS2 = as<NumericVector>(FITS2_SEXP);
  int EXITSTATUS = as<int>(EXITSTATUS_SEXP);

  // Call your C++ function
  BLACKBOXF(NRESPONDENTS, NISSUES, NDIMENSIONS, NMISSING,
            KMISS, MINSCALE, MID, KISSUE, FITS, PSIMATRIX,
            WMATRIX, LRESPONDENTS, LMARK, FITS2, EXITSTATUS);

  // Prepare the return value (if any)
  // If your function updates the passed variables, they are already updated in R
  // If you need to return a value, wrap it in a SEXP
  return R_NilValue; // or return appropriate SEXP if needed
}

// [[Rcpp::export]]
List blackboxt_wrapper(
    int NRESPONDENTS,
    int NISSUES,
    int NDIMENSIONS,
    int NMISSING,
    NumericMatrix KMISS,
    int MINSCALE,
    IntegerVector MID,
    NumericMatrix KISSUE,
    NumericVector FITS,
    NumericMatrix PSIMATRIX,
    NumericMatrix WMATRIX,
    IntegerVector LRESPONDENTS,
    IntegerVector LMARK,
    NumericVector FITS2,
    int EXITSTATUS
) {
  // Convert inputs to Eigen types if necessary
  Eigen::MatrixXd KMISS_Eigen = as<Eigen::MatrixXd>(KMISS);
  Eigen::VectorXi MID_Eigen = as<Eigen::VectorXi>(MID);
  Eigen::MatrixXd KISSUE_Eigen = as<Eigen::MatrixXd>(KISSUE);
  Eigen::VectorXd FITS_Eigen = as<Eigen::VectorXd>(FITS);
  Eigen::MatrixXd PSIMATRIX_Eigen = as<Eigen::MatrixXd>(PSIMATRIX);
  Eigen::MatrixXd WMATRIX_Eigen = as<Eigen::MatrixXd>(WMATRIX);
  Eigen::VectorXi LRESPONDENTS_Eigen = as<Eigen::VectorXi>(LRESPONDENTS);
  Eigen::VectorXi LMARK_Eigen = as<Eigen::VectorXi>(LMARK);
  Eigen::VectorXd FITS2_Eigen = as<Eigen::VectorXd>(FITS2);

  // Call your C++ function
  blackboxt(
    NRESPONDENTS,
    NISSUES,
    NDIMENSIONS,
    NMISSING,
    KMISS_Eigen,
    MINSCALE,
    MID_Eigen,
    KISSUE_Eigen,
    FITS_Eigen,
    PSIMATRIX_Eigen,
    WMATRIX_Eigen,
    LRESPONDENTS_Eigen,
    LMARK_Eigen,
    FITS2_Eigen,
    EXITSTATUS
  );

  // Prepare the return list
  return List::create(
    Named("FITS") = wrap(FITS_Eigen),
    Named("PSIMATRIX") = wrap(PSIMATRIX_Eigen),
    Named("WMATRIX") = wrap(WMATRIX_Eigen),
    Named("LRESPONDENTS") = wrap(LRESPONDENTS_Eigen),
    Named("LMARK") = wrap(LMARK_Eigen),
    Named("FITS2") = wrap(FITS2_Eigen),
    Named("EXITSTATUS") = EXITSTATUS
  );
}

// [[Rcpp::export]]
List mckalnew_wrapper(
    int NRESPONDENTS,
    int NISSUES,
    int NDIMENSIONS,
    int NSCALE,
    NumericMatrix KDATA,
    int MINSCALE,
    IntegerVector MID,
    NumericVector WEIGHTS,
    NumericVector VPARAMETERS,
    NumericVector KISSUE,
    NumericMatrix KMATRIX,
    double CONVERGE,
    int EXITSTATUS
) {
  // Convert inputs to Eigen types if necessary
  Eigen::MatrixXd KDATA_Eigen = as<Eigen::MatrixXd>(KDATA);
  Eigen::VectorXi MID_Eigen = as<Eigen::VectorXi>(MID);
  Eigen::VectorXd WEIGHTS_Eigen = as<Eigen::VectorXd>(WEIGHTS);
  Eigen::VectorXd VPARAMETERS_Eigen = as<Eigen::VectorXd>(VPARAMETERS);
  Eigen::VectorXd KISSUE_Eigen = as<Eigen::VectorXd>(KISSUE);
  Eigen::MatrixXd KMATRIX_Eigen = as<Eigen::MatrixXd>(KMATRIX);

  // Call your C++ function
  mckalnew(
    NRESPONDENTS,
    NISSUES,
    NDIMENSIONS,
    NSCALE,
    KDATA_Eigen,
    MINSCALE,
    MID_Eigen,
    WEIGHTS_Eigen,
    VPARAMETERS_Eigen,
    KISSUE_Eigen,
    KMATRIX_Eigen,
    CONVERGE,
    EXITSTATUS
  );

  // Prepare the return list
  return List::create(
    Named("WEIGHTS") = wrap(WEIGHTS_Eigen),
    Named("VPARAMETERS") = wrap(VPARAMETERS_Eigen),
    Named("KMATRIX") = wrap(KMATRIX_Eigen),
    Named("CONVERGE") = CONVERGE,
    Named("EXITSTATUS") = EXITSTATUS
  );
}
