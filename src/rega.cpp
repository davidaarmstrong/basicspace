#include <RcppEigen.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

// The REGA function using RcppEigen
// [[Rcpp::export]]
void REGA(int NISSUES, int NDIMENSIONS, int NS, int NF,
          const Eigen::MatrixXd& A, const Eigen::VectorXd& Y, Eigen::VectorXd& V, int IPRNT)
{
  // Suppress unused variable warning
  (void)IPRNT;

  // Allocate matrices and vectors
  Eigen::MatrixXd B(NF, NF);
  Eigen::MatrixXd ZMAT(NF, NF);
  Eigen::VectorXd WVEC(NF);

  // Compute B(J, JJ) = sum over I of A(I, J) * A(I, JJ)
  for (int J = 0; J < NF; ++J) {
    for (int JJ = 0; JJ < NF; ++JJ) {
      double SUM = 0.0;
      for (int I = 0; I < NS; ++I) {
        SUM += A(I, J) * A(I, JJ);
      }
      B(J, JJ) = SUM;
      ZMAT(J, JJ) = SUM;
    }
  }

  // Compute eigenvalues and eigenvectors of ZMAT
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ZMAT);
  if (es.info() != Eigen::Success) {
    Rcpp::Rcerr << "Eigenvalue decomposition failed." << std::endl;
    return;
  }
  WVEC = es.eigenvalues();
  Eigen::MatrixXd eigenvectors = es.eigenvectors();

  // Compute inverse of B using eigen decomposition
  Eigen::MatrixXd D_inv = Eigen::MatrixXd::Zero(NF, NF);
  for (int i = 0; i < NF; ++i) {
    if (std::abs(WVEC(i)) > 0.0001) {
      D_inv(i, i) = 1.0 / WVEC(i);
    }
  }
  Eigen::MatrixXd C = eigenvectors * D_inv * eigenvectors.transpose();

  // Compute BB(J, I) = sum over JJ of C(J, JJ) * A(I, JJ)
  Eigen::MatrixXd BB(NF, NS);
  for (int J = 0; J < NF; ++J) {
    for (int I = 0; I < NS; ++I) {
      double SUM = 0.0;
      for (int JJ = 0; JJ < NF; ++JJ) {
        SUM += C(J, JJ) * A(I, JJ);
      }
      BB(J, I) = SUM;
    }
  }

  // Compute V(JJ) = sum over J of BB(JJ, J) * Y(J)
  for (int JJ = 0; JJ < NF; ++JJ) {
    double SUM = 0.0;
    for (int J = 0; J < NS; ++J) {
      SUM += BB(JJ, J) * Y(J);
    }
    V(JJ) = SUM;
  }

  // No need to deallocate in C++, as Eigen automatically manages memory
}
