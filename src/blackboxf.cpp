#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// [[Rcpp::depends(RcppEigen)]]

// The BLACKBOXF function using RcppEigen
// [[Rcpp::export]]
void BLACKBOXF(
    int NRESPONDENTS,
    int NISSUES,
    int NDIMENSIONS,
    int NMISSING,
    const Eigen::MatrixXd& KMISS,          // KMISS(NISSUES, NMISSING)
    int MINSCALE,
    const Eigen::VectorXi& MID,            // MID(NRESPONDENTS)
    const Eigen::MatrixXd& KISSUE,         // KISSUE(NRESPONDENTS, NISSUES)
    Eigen::VectorXd& FITS,                 // FITS(7 * NDIMENSIONS)
    Eigen::MatrixXd& PSIMATRIX,            // PSIMATRIX (size depends on NDIMENSIONS)
    Eigen::MatrixXd& WMATRIX,              // WMATRIX (size depends on NDIMENSIONS)
    Eigen::VectorXi& LRESPONDENTS,         // LRESPONDENTS(NISSUES)
    Eigen::VectorXi& LMARK,                // LMARK(NRESPONDENTS)
    Eigen::VectorXd& FITS2,                // FITS2(6)
    int& EXITSTATUS
) {
  // Initialize variables
  int NS = NDIMENSIONS;
  int NY = NISSUES;
  int NMISS = NMISSING;
  int KVMIN = MINSCALE;
  int II = 0;
  int LWORK = 3 * NRESPONDENTS + 3 * NISSUES;
  int IPRNT = 1; // Printer switch: 1 - no write to disk, 0 - write to disk

  // Allocate matrices and vectors
  Eigen::VectorXi KID(NRESPONDENTS);
  Eigen::VectorXi LID(NISSUES);
  Eigen::VectorXd YHAT(NRESPONDENTS);
  Eigen::VectorXd WORK(LWORK);
  Eigen::MatrixXd XBIGONE(NRESPONDENTS, NISSUES);
  Eigen::MatrixXd W(NISSUES, NDIMENSIONS + 2);
  Eigen::MatrixXd XDATA(NRESPONDENTS, NDIMENSIONS);
  Eigen::MatrixXd XT(NRESPONDENTS, NISSUES);
  Eigen::MatrixXd RSAVE(NDIMENSIONS, 3);
  Eigen::MatrixXd UUU(NRESPONDENTS, NISSUES);
  Eigen::MatrixXd VVV(NRESPONDENTS, NISSUES);

  // Initialize LMARK
  LMARK.setZero(NRESPONDENTS);

  // Process data and handle missing values
  for (int I = 0; I < NRESPONDENTS; ++I) {
    int IMARK = 0;
    for (int J = 0; J < NY; ++J) {
      for (int K = 0; K < NMISS; ++K) {
        if (std::abs(KISSUE(I, J) - KMISS(J, K)) <= 0.001) {
          KISSUE(I, J) = -999.0;
          IMARK++;
          break;
        }
      }
    }
    if ((NY - IMARK) < KVMIN)
      continue;
    LMARK[I] = 1;
    KID[II] = MID[I];
    XBIGONE.row(II) = KISSUE.row(I);
    II++;
  }
  int NOBS = II;

  // Variables to keep track of positions in output matrices
  int KKKK = 0;
  int KKKKK = 0;

  // Loop over dimensions
  for (int KKK = 1; KKK <= NS; ++KKK) {
    // Call the BLACKB subroutine (assumed to be implemented elsewhere)
    double SVSUM = 0.0;
    FITS2.resize(6);
    // Assuming BLACKB is another function that you have defined elsewhere
    BLACKB(NOBS, NRESPONDENTS, NISSUES, NDIMENSIONS, KKK,
           XBIGONE.topRows(NOBS), XDATA.topRows(NOBS), W, SVSUM, FITS2, IPRNT);

    // Constraint checks and computations
    // Compute sums for constraint checks
    for (int K = 0; K < KKK; ++K) {
      double SUM = XDATA.col(K).sum();
      VVV(0, K) = SUM / static_cast<double>(NOBS);
    }

    // Compute cross-products of XDATA
    for (int J = 0; J < KKK; ++J) {
      for (int K = 0; K < KKK; ++K) {
        double SUM = (XDATA.col(J).array() * XDATA.col(K).array()).sum();
        VVV(J, K) = SUM;
      }
    }

    // Compute R-squares and other statistics
    double ASUM = 0.0, BSUM = 0.0, CSUM = 0.0, DSUM = 0.0, ESUM = 0.0, SUME = 0.0;
    int KK = 0;
    for (int J = 0; J < NY; ++J) {
      double AASUM = 0.0, BBSUM = 0.0, CCSUM = 0.0, DDSUM = 0.0, EESUM = 0.0;
      int KJJ = 0;
      for (int I = 0; I < NOBS; ++I) {
        double SUM = 0.0;
        for (int K = 0; K < KKK; ++K) {
          SUM += XDATA(I, K) * W(J, K + 1);
        }
        XT(I, J) = SUM;
        double AA = SUM + W(J, 0); // W(J,1) in Fortran is W(J,0) in C++
        if (std::abs(XBIGONE(I, J) + 999.0) <= 0.001)
          continue;
        double BB = XBIGONE(I, J);
        SUME += std::pow(AA - BB, 2);
        AASUM += AA;
        BBSUM += BB;
        CCSUM += AA * AA;
        DDSUM += BB * BB;
        EESUM += AA * BB;
        KJJ++;
      }
      double AAA = static_cast<double>(KJJ) * EESUM - AASUM * BBSUM;
      double BBB = static_cast<double>(KJJ) * CCSUM - AASUM * AASUM;
      double CCC = static_cast<double>(KJJ) * DDSUM - BBSUM * BBSUM;
      double RRR = 0.0;
      if (std::abs(BBB * CCC) > 0.0)
        RRR = (AAA * AAA) / (BBB * CCC);
      W(J, KKK + 1) = RRR; // Store R-square
      LID[J] = KJJ;
      ASUM += AASUM;
      BSUM += BBSUM;
      CSUM += CCSUM;
      DSUM += DDSUM;
      ESUM += EESUM;
      KK += KJJ;
    }
    double AA = static_cast<double>(KK) * ESUM - ASUM * BSUM;
    double BB = static_cast<double>(KK) * CSUM - ASUM * ASUM;
    double CC = static_cast<double>(KK) * DSUM - BSUM * BSUM;
    double RRR = (AA * AA) / (BB * CC);
    RSAVE(KKK - 1, 0) = SUME;
    RSAVE(KKK - 1, 1) = RRR;
    RSAVE(KKK - 1, 2) = std::sqrt(SUME / static_cast<double>(KK - KKK * (NOBS + NY) - NY));

    // Singular Value Decomposition (SVD)
    Eigen::JacobiSVD<MatrixXd> svd(XT.topRows(NOBS), Eigen::ComputeThinU | Eigen::ComputeThinV);
    YHAT.head(NOBS) = svd.singularValues();

    // Write out column parameters and transfer to WMATRIX
    for (int I = 0; I < NY; ++I) {
      for (int JJ = 0; JJ < KKK + 2; ++JJ) {
        if (KKKK == 0) {
          WMATRIX(I * (KKK + 2) + JJ) = W(I, JJ);
        } else {
          WMATRIX((KKKK) * NY + I * (KKK + 2) + JJ) = W(I, JJ);
        }
      }
      LRESPONDENTS[I] = LID[I];
    }

    // Write out row parameters and transfer to PSIMATRIX
    int II = 0;
    for (int I = 0; I < NRESPONDENTS; ++I) {
      if (LMARK[I] == 0)
        continue;
      for (int JJ = 0; JJ < KKK; ++JJ) {
        if (KKKKK == 0) {
          PSIMATRIX(I * KKK + JJ) = XDATA(II, JJ);
        } else {
          PSIMATRIX(KKKKK * NRESPONDENTS + I * KKK + JJ) = XDATA(II, JJ);
        }
      }
      II++;
    }

    // Update positions for output matrices
    KKKK = KKKK + KKK + 2;
    KKKKK = KKKKK + KKK;
  }

  // Update FITS array
  for (int J = 0; J < NS; ++J) {
    double AA = (SVSUM - RSAVE(J, 0)) / SVSUM;
    double CCC = SVSUM - RSAVE(J, 0);
    double BB = (J == 0) ? SVSUM : RSAVE(J - 1, 0);
    double SUM = (BB - RSAVE(J, 0)) / SVSUM;

    FITS(0 + J * 7) = RSAVE(J, 0);
    FITS(1 + J * 7) = CCC;
    FITS(2 + J * 7) = SUM * 100.0;
    FITS(3 + J * 7) = AA * 100.0;
    FITS(4 + J * 7) = RSAVE(J, 1);
    FITS(5 + J * 7) = RSAVE(J, 2);
    FITS(6 + J * 7) = YHAT(J);
  }

  // Set exit status to success
  EXITSTATUS = 1;
}
