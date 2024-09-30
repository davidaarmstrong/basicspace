#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;  // For dynamic matrices
using Eigen::VectorXd;  // For dynamic vectors

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
void REG(int NRESPONDENTS, int NISSUES, int NDIMENSIONS, int NP, int NF, int NY,
         std::vector<double>& TSUM,
         Eigen::MatrixXd& W,
         Eigen::MatrixXd& XS,
         Eigen::MatrixXd& X,
         Eigen::MatrixXd& PSI,
         int IPRNT, int ILAST, int KKK, double& AREG, double& BREG)
{
  // Declare and allocate variables
  std::vector<int> LLL(NISSUES, 0);
  std::vector<int> LL(NISSUES, 0);
  std::vector<int> MPOS(NISSUES, 0);
  std::vector<double> C(NISSUES, 0.0);
  std::vector<double> RSUM(NISSUES, 0.0);

  Eigen::VectorXd WVEC(3 * NDIMENSIONS);
  Eigen::VectorXd WK(40 * NDIMENSIONS);
  Eigen::MatrixXd A(3 * NDIMENSIONS, 3 * NDIMENSIONS);
  Eigen::MatrixXd B(3 * NDIMENSIONS, 3 * NDIMENSIONS);
  Eigen::MatrixXd XT(NRESPONDENTS, NISSUES);
  Eigen::MatrixXd R(NISSUES, NISSUES);
  Eigen::MatrixXd ZMAT(3 * NDIMENSIONS, 3 * NDIMENSIONS);

  double H1947 = BREG;
  int KPX1947 = KKK;

  int NF1 = NF + 1;
  int LWORK = 40 * NDIMENSIONS;

  // Initialize variables
  for (int K = 0; K < NY; ++K) {
    LLL[K] = 0;
    TSUM[K] = 0.0;
    C[K] = 0.0;
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
      Eigen::MatrixXd PP = Eigen::MatrixXd::Zero(NF1, NF1);
      for (int J = 0; J < NF1; ++J) {
        for (int JJ = 0; JJ < NF1; ++JJ) {
          double SUM = 0.0;
          for (int I = 0; I < NP; ++I) {
            if (std::abs(XS(I, K) + 999.0) <= 0.001)
              continue;
            SUM += PSI(I, J) * PSI(I, JJ);
          }
          A(J, JJ) = SUM;
          ZMAT(J, JJ) = SUM;
        }
      }

      // Compute eigenvalues and eigenvectors of ZMAT
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(ZMAT.topLeftCorner(NF1, NF1));
      if (es.info() != Eigen::Success) {
        Rcpp::Rcerr << "Eigenvalue decomposition failed." << std::endl;
        return;
      }
      WVEC.head(NF1) = es.eigenvalues();
      MatrixXd eigenvectors = es.eigenvectors();

      // Compute inverse of A using eigen decomposition
      MatrixXd D_inv = MatrixXd::Zero(NF1, NF1);
      for (int i = 0; i < NF1; ++i) {
        if (std::abs(WVEC(i)) > 0.0001) {
          D_inv(i, i) = 1.0 / WVEC(i);
        }
      }
      B.topLeftCorner(NF1, NF1) = eigenvectors * D_inv * eigenvectors.transpose();

      // Compute [P'P]^{-1} * P' * X = W'/c'
      for (int I = 0; I < NP; ++I) {
        if (std::abs(XS(I, K) + 999.0) <= 0.001)
          continue;

        // Compute C vector
        Eigen::VectorXd C_vec = B.topLeftCorner(NF1, NF1) * PSI.row(I).head(NF1).transpose();

        // Update W matrix
        for (int J = 0; J < NF1; ++J) {
          W(K, J) += C_vec(J) * XS(I, K);
        }
      }

      // Calculate R-square and SSE for regression
      double ESUM = 0.0;
      for (int I = 0; I < NP; ++I) {
        if (std::abs(XS(I, K) + 999.0) <= 0.001)
          continue;
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

      int NYY = 2;
      // Initialize LL and MPOS vectors
      std::fill(LL.begin(), LL.end(), 0);
      std::fill(MPOS.begin(), MPOS.end(), 0);
      int KS = 0, KPOS = 0;

      // Call external function CORR22 (Placeholder: implement the function or add as needed)
      CORR22(NRESPONDENTS, NISSUES, NP, NYY, XT.leftCols(2), R, LL, MPOS, KS, KPOS, 1);

      RSUM[K] = R(0, 1) * R(0, 1);
    }

    // End of regression loop for current dimension
    if (ILAST == 1) {
      goto end_loop;
    }

    // If ILAST = 0, calculate the R-square for the current number of dimensions
    for (int I = 0; I < NP; ++I) {
      for (int K = 0; K < NY; ++K) {
        double SUM = 0.0;
        for (int J = 0; J < NF1; ++J) {
          SUM += PSI(I, J) * W(K, J);
        }
        XT(I, K) = SUM;
      }
    }
    double RR = 0.0;

    // Call external function RSQUR (Placeholder: implement the function or add as needed)
    RSQUR(NRESPONDENTS, NISSUES, NP, NY, RR, XT, XS, 1);
    TSUM[KK - 1 + NY] = RR;

    end_loop:
      if (IPRNT == 1)
        continue;

      if (ILAST == 1) {
        AREG = TSUM[KK - 1];
      } else if (ILAST == 0) {
        AREG = TSUM[KK - 1];
      }
  }

  // Variables are automatically deallocated when going out of scope
}
