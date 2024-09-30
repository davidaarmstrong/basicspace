#include <RcppEigen.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

// The REGT function using RcppEigen
// [[Rcpp::export]]
void REGT(int NRESPONDENTS, int NISSUES, int NDIMENSIONS, int NP, int NF, int NY,
          std::vector<double>& TSUM,    // TSUM(2*NRESPONDENTS)
          Eigen::MatrixXd& W,           // W(NRESPONDENTS, NDIMENSIONS+2)
          Eigen::MatrixXd& XS,          // XS(NISSUES, NRESPONDENTS)
          Eigen::MatrixXd& X,           // X(NISSUES, NRESPONDENTS)
          Eigen::MatrixXd& PSI,         // PSI(NISSUES, NRESPONDENTS)
          int IPRNT, int ILAST, int KKK, double& AREG)
{
  // Suppress unused variable warnings
  (void)IPRNT;
  (void)KKK;

  int NF1 = NF + 1;
  int LWORK = 40 * (NDIMENSIONS + 2);

  // Allocate vectors and matrices
  std::vector<int> LLL(NY, 0);
  VectorXd C(NY);
  VectorXd WVEC(NF1);
  MatrixXd A(NF1, NF1);
  MatrixXd B(NF1, NF1);
  MatrixXd XT(NP, NY); // Transposed matrices compared to Fortran code
  MatrixXd R(NY, NY);
  MatrixXd ZMAT(NF1, NF1);

  // Initialize variables
  for (int K = 0; K < NY; ++K) {
    LLL[K] = 0;
    TSUM[K] = 0.0;
    C(K) = 0.0;
    for (int J = 0; J < NF1; ++J) {
      W(K, J) = 0.0;
    }
  }

  // Determine starting dimension based on ILAST
  int KKA = (ILAST == 1) ? NF : 1;

  // Loop over dimensions
  for (int KK = KKA; KK <= NF; ++KK) {
    // Initialize W matrix
    for (int K = 0; K < NY; ++K) {
      for (int J = 0; J <= NF; ++J) {
        W(K, J) = 0.0;
      }
    }

    NF1 = KK + 1;
    int NF11 = KK + 2;

    // Loop over issues
    for (int K = 0; K < NY; ++K) {
      LLL[K] = 0;

      // Compute [P'P] matrix
      MatrixXd PP = MatrixXd::Zero(NF1, NF1);
      for (int J = 0; J < NF1; ++J) {
        for (int JJ = 0; JJ < NF1; ++JJ) {
          double SUM = 0.0;
          for (int I = 0; I < NP; ++I) {
            // Check for missing data
            if (std::abs(XS(I, K) + 999.0) <= 0.001) continue;
            if (std::abs(PSI(I, 0) + 999.0) <= 0.001) continue;
            SUM += PSI(I, J) * PSI(I, JJ);
          }
          A(J, JJ) = SUM;
          ZMAT(J, JJ) = SUM;
        }
      }

      // Compute eigenvalues and eigenvectors of ZMAT
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(ZMAT);
      if (es.info() != Eigen::Success) {
        Rcpp::Rcerr << "Eigenvalue decomposition failed." << std::endl;
        return;
      }
      WVEC = es.eigenvalues();
      MatrixXd eigenvectors = es.eigenvectors();

      // Compute inverse of A using eigen decomposition
      MatrixXd D_inv = MatrixXd::Zero(NF1, NF1);
      for (int i = 0; i < NF1; ++i) {
        if (std::abs(WVEC(i)) > 0.001) {
          D_inv(i, i) = 1.0 / WVEC(i);
        }
      }
      B = eigenvectors * D_inv * eigenvectors.transpose();

      // Compute [P'P]^{-1} * P' * X = W'/c'
      for (int I = 0; I < NP; ++I) {
        // Check for missing data
        if (std::abs(XS(I, K) + 999.0) <= 0.001) continue;
        if (std::abs(PSI(I, 0) + 999.0) <= 0.001) continue;

        // Compute C vector
        VectorXd PSI_row = PSI.row(I).head(NF1);
        VectorXd C_vec = B * PSI_row;

        // Update W matrix
        for (int J = 0; J < NF1; ++J) {
          W(K, J) += C_vec(J) * XS(I, K);
        }
      }

      // Calculate R-square and SSE for regression
      double ESUM = 0.0;
      XT.col(0).setConstant(-999.0);
      XT.col(1).setConstant(-999.0);
      for (int I = 0; I < NP; ++I) {
        // Check for missing data
        if (std::abs(XS(I, K) + 999.0) <= 0.001) continue;
        if (std::abs(PSI(I, 0) + 999.0) <= 0.001) continue;
        LLL[K] += 1;
        double SUM = 0.0;
        for (int J = 0; J < NF1; ++J) {
          SUM += PSI(I, J) * W(K, J);
        }
        ESUM += std::pow(SUM - XS(I, K), 2);
        X(I, K) = SUM - XS(I, K);
        XT(I, 0) = SUM;
        XT(I, 1) = XS(I, K);
      }
      W(K, NF11 - 1) = ESUM;
      TSUM[KK - 1] += ESUM;
    }

    // End of regression loop for current dimension
    if (ILAST == 1) {
      AREG = TSUM[KK - 1];
      continue;
    }

    // If ILAST == 0, calculate the R-square for the current number of dimensions
    XT.setConstant(-999.0);
    for (int I = 0; I < NP; ++I) {
      // Check for missing data
      if (std::abs(PSI(I, 0) + 999.0) <= 0.001) continue;
      for (int K = 0; K < NY; ++K) {
        double SUM = 0.0;
        for (int J = 0; J < NF1; ++J) {
          SUM += PSI(I, J) * W(K, J);
        }
        XT(I, K) = SUM;
      }
    }

    // Call RSQUR function
    double RR = 0.0;
    RSQUR(NISSUES, NRESPONDENTS, NP, NY, RR, XT, XS, IPRNT);
    TSUM[KK - 1 + NY] = RR;

    AREG = TSUM[KK - 1];
  }

  // Variables are automatically deallocated when going out of scope
}
