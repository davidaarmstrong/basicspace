#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;  // For dynamic matrices

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
void RSQUR(
    int NRESPONDENTS,
    int NISSUES,
    int NP,
    int NY,
    double& R,
    const Eigen::MatrixXd& A,  // Use Eigen::MatrixXd for A
    const Eigen::MatrixXd& B,  // Use Eigen::MatrixXd for B
    int IPRNT)
{
  // Variables for accumulating sums
  double ASUM = 0.0;
  double BSUM = 0.0;
  double CSUM = 0.0;
  double DSUM = 0.0;
  double ESUM = 0.0;
  double XNT = 0.0;

  // Loop over each data point
  for (int I = 0; I < NP; ++I) {
    for (int J = 0; J < NY; ++J) {
      // Skip if A(I, J) or B(I, J) is approximately equal to -999.0
      if (std::abs(A(I, J) + 999.0) <= 0.001) continue;
      if (std::abs(B(I, J) + 999.0) <= 0.001) continue;

      // Accumulate sums
      ASUM += A(I, J);
      BSUM += B(I, J);
      CSUM += A(I, J) * A(I, J);
      DSUM += B(I, J) * B(I, J);
      ESUM += A(I, J) * B(I, J);
      XNT += 1.0;
    }
  }

  // Compute intermediate values
  double AA = XNT * ESUM - ASUM * BSUM;
  double BB = XNT * CSUM - ASUM * ASUM;
  double CC = XNT * DSUM - BSUM * BSUM;

  // Compute the squared correlation coefficient R
  if (BB * CC != 0.0) {
    R = (AA * AA) / (BB * CC);
  } else {
    // Handle division by zero if necessary
    R = 0.0;
  }
}
