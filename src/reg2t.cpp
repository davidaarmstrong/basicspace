#include <RcppEigen.h>
#include <iostream>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Function prototype for REGAT (to be implemented elsewhere)
void REGAT(int NRESPONDENTS, int NDIMENSIONS, int NS, int NF, const Eigen::MatrixXd& A_sub, const Eigen::VectorXd& Y_sub, Eigen::VectorXd& V_sub);

// The REG2T function using RcppEigen
// [[Rcpp::export]]
void REG2T(int NRESPONDENTS, int NISSUES, int NDIMENSIONS, int NP, int NF, int NY,
           Eigen::MatrixXd& W,     // W(NISSUES, NDIMENSIONS+2)
           Eigen::MatrixXd& XS,    // XS(NISSUES, NRESPONDENTS)
           Eigen::MatrixXd& X,     // X(NISSUES, NRESPONDENTS)
           Eigen::MatrixXd& PSI,   // PSI(NISSUES, NRESPONDENTS)
           double& PXB,
           double& PXS,
           int KKK,
           int NWHO,
           double& BREG)
{
  int NF1 = NWHO + 1;
  double ESUM = 0.0;
  PXB = 0.0;
  PXS = 0.0;
  double XNS = 0.0;

  // Allocate temporary variables
  Eigen::VectorXd V_vec(NF);
  Eigen::VectorXd Y_vec(NRESPONDENTS);            // Maximum possible size
  Eigen::MatrixXd A_mat(NRESPONDENTS, NF);        // Maximum possible size

  // Regress one respondent at a time because of missing data
  for (int I = 0; I < NP; ++I)
  {
    // Check if PSI(I, 0) is missing data
    if (std::abs(PSI(I, 0) + 999.0) <= 0.001)
      continue;

    int KK = 0;

    // Set up dependent variable Y and matrix A
    for (int J = 0; J < NY; ++J)
    {
      // Check for missing data in XS
      if (std::abs(XS(J, I) + 999.0) <= 0.001)
        continue;

      Y_vec(KK) = XS(J, I) - W(J, NF1 - 1);  // Adjusted indices
      for (int JJ = 0; JJ < NF; ++JJ)
      {
        A_mat(KK, JJ) = W(J, JJ);
      }
      KK++;
    }

    int NS = KK;

    if (NS == 0)
    {
      // No valid data, skip to next respondent
      continue;
    }

    // Resize A_mat and Y_vec to actual size
    Eigen::MatrixXd A_sub = A_mat.topRows(NS);
    Eigen::VectorXd Y_sub = Y_vec.head(NS);

    // Call REGAT to perform least squares
    Eigen::VectorXd V_sub(NF);
    REGAT(NRESPONDENTS, NDIMENSIONS, NS, NF, A_sub, Y_sub, V_sub);

    // Store estimated entries of PSI and compute SSE
    for (int K = 0; K < NY; ++K)
    {
      double SUM = 0.0;
      if (std::abs(XS(K, I) + 999.0) <= 0.001)
        continue;

      for (int J = 0; J < NF; ++J)
      {
        PSI(K, I) = V_sub(J);
        SUM += PSI(K, I) * W(K, J);
      }
      SUM += W(K, NF1 - 1); // W(K, NF1)

      X(K, I) = SUM - XS(K, I);
      ESUM += std::pow(SUM - XS(K, I), 2);
    }

    PXB += PSI(0, I);
    PXS += std::pow(PSI(0, I), 2);
    XNS += 1.0;
  }

  // Compute PXB and PXS
  if (XNS > 0)
  {
    PXB = PXB / XNS;
    PXS = PXS - XNS * PXB * PXB;
  }
  else
  {
    PXB = 0.0;
    PXS = 0.0;
  }

  BREG = ESUM;
}
