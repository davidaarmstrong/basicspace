#include <RcppEigen.h>
#include <iostream>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// The REG2 function using RcppEigen
// [[Rcpp::export]]
void REG2(int NRESPONDENTS, int NISSUES, int NDIMENSIONS, int NP, int NF, int NY,
          Eigen::MatrixXd& W,  // W(NISSUES, NDIMENSIONS+2)
          Eigen::MatrixXd& XS, // XS(NRESPONDENTS, NISSUES)
          Eigen::MatrixXd& X,  // X(NRESPONDENTS, NISSUES)
          Eigen::MatrixXd& PSI,// PSI(NRESPONDENTS, NISSUES)
          double& PXB, double& PXS,
          int KKK, int IPRNT, double AREG, double& BREG)
{
  // Suppress unused variable warnings
  (void)IPRNT;
  (void)KKK;
  (void)AREG;

  int NF1 = NF + 1;
  double ESUM = 0.0;
  PXB = 0.0;
  PXS = 0.0;
  double XNS = 0.0;

  // Allocate temporary variables
  VectorXd V(NF);
  VectorXd Y(NISSUES);
  MatrixXd A_matrix(NISSUES, NF);

  // Regress one row at a time due to missing data
  for (int I = 0; I < NP; ++I) {
    int KK = 0;

    // Set up dependent variable Y and matrix A
    for (int J = 0; J < NY; ++J) {
      V(J) = 0.0;

      // Check for missing data (value approximately equal to -999.0)
      if (std::abs(XS(I, J) + 999.0) <= 0.001)
        continue;

      // Increment valid data count
      Y(KK) = XS(I, J) - W(J, NF1 - 1); // W indices adjusted for zero-based indexing
      for (int JJ = 0; JJ < NF; ++JJ) {
        A_matrix(KK, JJ) = W(J, JJ);
      }
      KK++;
    }

    int NS = KK;

    if (NS == 0) {
      // No valid data for this respondent, skip to next
      continue;
    }

    // Resize A and Y to match the number of valid data points
    MatrixXd A_sub = A_matrix.topRows(NS);
    VectorXd Y_sub = Y.head(NS);

    // Perform least squares estimation using Eigen's linear solver
    VectorXd V_sub = A_sub.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y_sub);

    // Store estimated row entries of PSI and compute SSE
    for (int K = 0; K < NY; ++K) {
      double SUM = 0.0;

      // Check for missing data
      if (std::abs(XS(I, K) + 999.0) <= 0.001)
        continue;

      for (int J = 0; J < NF; ++J) {
        PSI(I, J) = V_sub(J);
        SUM += PSI(I, J) * W(K, J);
      }

      SUM += W(K, NF1 - 1); // Adjusted indexing for zero-based

      X(I, K) = SUM - XS(I, K);
      ESUM += std::pow(SUM - XS(I, K), 2);
    }

    PXB += PSI(I, 0);
    PXS += std::pow(PSI(I, 0), 2);
    XNS += 1.0;
  }

  // Compute PXB and PXS
  if (XNS != 0.0) {
    PXB = PXB / XNS;
    PXS = PXS - XNS * PXB * PXB;
  } else {
    PXB = 0.0;
    PXS = 0.0;
  }

  BREG = ESUM;

  // No need to deallocate in C++, as vectors and matrices are automatically managed
}
