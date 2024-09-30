#include <RcppEigen.h>
#include <iostream>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;

// The MCKALNEW function using RcppEigen
// [[Rcpp::export]]
void MCKALNEW(
    int NRESPONDENTS,
    int NISSUES,
    int NSELFPOS,
    int NMISSING,
    const Eigen::VectorXd& KMISS,           // KMISS(NISSUES * NMISSING)
    int POLARITY,
    const Eigen::VectorXi& MID,             // MID(NRESPONDENTS)
    const Eigen::MatrixXd& KISSUE,          // KISSUE(NRESPONDENTS, NISSUES + 1)
    Eigen::VectorXd& FITS,                  // FITS(5)
    Eigen::MatrixXd& PSIMATRIX,             // PSIMATRIX(NRESPONDENTS, 4)
    Eigen::VectorXd& STIMCOORDS,            // STIMCOORDS(NISSUES)
    Eigen::VectorXd& EIGENVALUES,           // EIGENVALUES(NISSUES)
    int& EXITSTATUS
) {
  // Declare and allocate variables
  int NQ = NISSUES - 1;
  int NRESP = NSELFPOS;
  if (NRESP == 0)
    NQ = NISSUES;

  int NNQ = NQ + 1;
  if (NRESP == 0)
    NNQ = NQ;

  int NLEFT = POLARITY;
  if (NLEFT == 0)
    NLEFT = 1;

  // Initialize arrays and matrices
  Eigen::MatrixXi KDATA(NRESPONDENTS, 2); // KDATA(NRESPONDENTS, 2)
  Eigen::VectorXd KDEGO(NQ);              // KDEGO(3 * NISSUES)
  Eigen::VectorXd A(NISSUES * NISSUES);
  Eigen::VectorXd ZZ(3 * NISSUES);
  Eigen::VectorXd XMISS(3 * NMISSING);
  Eigen::VectorXd B(3 * NISSUES);
  Eigen::VectorXd D(6 * NISSUES);
  Eigen::VectorXd WK(6 * NISSUES);
  Eigen::VectorXd ZSUM(3 * NISSUES);
  MatrixXd AWV(NQ, NQ);                   // AWV(NQ, NQ)
  MatrixXd XDATA(NRESPONDENTS, NNQ);      // XDATA(NRESPONDENTS, NNQ)

  // Initialize variables
  int KERR = 0;
  double XSIGM = 0.0;
  double BTSUM = 0.0;
  double ZZSUM = 0.0;
  int NP = 0;
  double XNQ = static_cast<double>(NQ);
  double XKNTL = 0.01;
  int KDEGN = 0;
  int KDEGM = 0;

  // Handle missing data
  if (NMISSING > 0) {
    for (int JJ = 0; JJ < 3 * NMISSING; ++JJ) {
      XMISS[JJ] = KMISS[JJ];
    }
  }

  // Initialize ZSUM and AWV matrices
  ZSUM.setZero();
  AWV.setZero();

  // Begin main loop
  for (int I = 0; I < NRESPONDENTS; ++I) {
    // Read data from KISSUE and MID
    Eigen::VectorXd A_row = KISSUE.row(I).head(NNQ);
    int I1 = MID[I];

    // Check for degeneracy
    double SUM = 0.0;
    int KXX = (NRESP == 1) ? 2 : 1;
    for (int JJ = 0; JJ < NNQ; ++JJ) {
      B[JJ] = A_row[JJ];
      if (NRESP == JJ + 1)
        continue;
      SUM += std::abs(A_row[JJ] - A_row[KXX - 1]);
    }
    if (SUM > XKNTL) {
      goto PROCESS_DATA;
    }
    if (NMISSING == 0) {
      KDEGN++;
      continue;
    } else {
      for (int II = 0; II < NMISSING; ++II) {
        if (A_row[KXX - 1] == XMISS[II]) {
          KDEGM++;
          goto NEXT_RESPONDENT;
        }
      }
      KDEGN++;
      continue;
    }

    PROCESS_DATA:
      if (NMISSING != 0) {
        KERR = 0;
        for (int II = 0; II < NMISSING; ++II) {
          for (int K = 0; K < NNQ; ++K) {
            if (K == NRESP - 1)
              continue;
            double AA = std::abs(A_row[K] - XMISS[II]);
            if (AA <= XKNTL) {
              KDEGO[K]++;
              KERR = 1;
            }
          }
        }
        if (KERR == 1)
          continue;
      }

      // Store data
      KDATA(NP, 0) = I;
      KDATA(NP, 1) = I1;
      XDATA.row(NP) = A_row.transpose();
      NP++;

      // Process data
      double ASUM = 0.0;
      double BSUM = 0.0;
      for (int II = 0; II < NQ; ++II) {
        ASUM += B[II];
        BSUM += B[II] * B[II];
        ZSUM[II] += B[II];
      }
      ZZSUM += ASUM;

      double DET = 1.0 / (XNQ * BSUM - ASUM * ASUM);
      for (int II = 0; II < NQ; ++II) {
        for (int J = 0; J <= II; ++J) {
          AWV(II, J) += DET * (BSUM - B[II] * ASUM + XNQ * B[II] * B[J] - B[J] * ASUM);
        }
      }

      NEXT_RESPONDENT:
        continue;
  }

  double XNP = static_cast<double>(NP);

  // Compute statistics
  ZZSUM = ZZSUM / (XNQ * XNP);
  for (int I = 0; I < NQ; ++I) {
    AWV(I, I) -= XNP;
    ZSUM[I] = (ZSUM[I] / XNP) - ZZSUM;
    BTSUM += ZSUM[I] * ZSUM[I];
  }
  BTSUM = std::sqrt(BTSUM);

  // Normalize ZSUM
  ZSUM /= BTSUM;

  // Perform eigenvalue decomposition
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(AWV);
  if (es.info() != Eigen::Success) {
    Rcpp::Rcout << "Eigenvalue decomposition failed." << std::endl;
    EXITSTATUS = -1;
    return;
  }
  Eigen::VectorXd D = es.eigenvalues();
  Eigen::MatrixXd eigenvectors = es.eigenvectors();

  // Transfer eigenvalues
  EIGENVALUES = D.head(NQ);

  double XROOT = D(NQ - 2); // Adjusted index for zero-based indexing
  double A1 = XNP * XROOT * (-1.0);
  double A2 = XNQ * std::pow(XNP + XROOT, 2);
  XSIGM = A1 / A2;

  // Transfer corrected goodness of fit
  FITS[0] = XSIGM;

  // Adjust eigenvectors based on POLARITY
  if (eigenvectors(NLEFT - 1, NQ - 1) > 0.0) {
    eigenvectors.col(NQ - 1) = -eigenvectors.col(NQ - 1);
  }

  // Transfer stimulus coordinates
  STIMCOORDS = eigenvectors.col(NQ - 1) * (XNP / (XNP + XROOT));

  // Call AMREG function (to be implemented)
  AMREG(NRESPONDENTS, NP, NQ, NRESP, NMISSING, ZZ, XMISS, BTSUM, ZZSUM, ZSUM, XSIGM, KDATA, XDATA, FITS, PSIMATRIX);

  // Set EXITSTATUS
  EXITSTATUS = 1;
}
