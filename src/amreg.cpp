#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;  // Explicitly use Eigen's MatrixXd
using Eigen::VectorXd;  // Explicitly use Eigen's VectorXd

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
void AMREG(
    int NRESPONDENTS,
    int NN,
    int NQ,
    int NRESP,
    int NMISS,
    const std::vector<double>& ZZ,         // ZZ(3*NQ)
    const std::vector<double>& XMISS,      // XMISS(3*NMISS)
    double BTSUM,
    double ZZSUM,
    const std::vector<double>& ZTSUM,      // ZTSUM(3*NQ)
    double XSIGM,
    const Eigen::MatrixXi& KDATA,          // KDATA(NRESPONDENTS,2)
    const Eigen::MatrixXd& XDATA,          // XDATA(NRESPONDENTS,NQ+1)
    std::vector<double>& FITS,             // FITS(5)
    std::vector<double>& PSIMATRIX         // PSIMATRIX(NRESPONDENTS*4)
)
{
  // Allocate arrays
  int NNQ = NQ + 1;
  int XNQ = NQ;

  // Vectors and matrices
  std::vector<double> A(NNQ, 0.0);
  std::vector<double> XSUM(2 * NQ, 0.0);
  std::vector<double> YSUM(NNQ, 0.0);
  std::vector<double> ZSUM(NNQ, 0.0);
  std::vector<double> Q(NN + NQ, 0.0);
  std::vector<double> QQ(NN + NQ, 0.0);
  std::vector<double> QQQ(NN + NQ, 0.0);
  MatrixXd XV(NQ, 2);    // Now explicitly Eigen::MatrixXd
  MatrixXd YV(2, 2);     // Now explicitly Eigen::MatrixXd
  MatrixXd ZV(2, NQ);    // Now explicitly Eigen::MatrixXd

  // Initialize variables
  double XSIGM2 = XSIGM;
  int NP = 0;
  double XKNTL = 0.001;
  double XMAX = -99.0;
  double TMAX = -99.0;
  double ASUM = 0.0;
  double BSUM = 0.0;
  double CSUM = 0.0;
  double DSUM = 0.0;
  double ESUM = 0.0;
  double TSUM = 0.0;
  double TTSUM = 0.0;

  // Initialize arrays
  for (int I = 0; I < NQ; ++I) {
    TMAX = std::max(TMAX, std::abs(ZZ[I]));
    Q[I] = ZZ[I];
    QQ[I] = ZZ[I];
    QQQ[I] = ZZ[I];
    XSUM[I] = 0.0;
    ZSUM[I] = 0.0;
    XSUM[I + NQ] = 0.0;
  }

  int IK = 0;
  int IPOS = 0;
  int INEG = 0;
  double SSE = 0.0;

  // Main loop
  while (true) {
    IK++;
    if (IK > NN)
      break;

    // Read data from KDATA and XDATA
    int IC = KDATA(IK - 1, 0);
    int I1 = KDATA(IK - 1, 1);
    for (int J = 0; J < NNQ; ++J) {
      A[J] = XDATA(IK - 1, J);
    }

    double XRES1 = 0.0;
    if (NRESP != 0) {
      if (NMISS != 0) {
        bool skip = false;
        for (int I = 0; I < NMISS; ++I) {
          double AA = std::abs(A[NRESP - 1] - XMISS[I]);
          if (AA <= XKNTL) {
            skip = true;
            break;
          }
        }
        if (skip)
          continue;
      }
      XRES1 = A[NRESP - 1];
    }

    if (NRESP == 0)
      XRES1 = 0.0;

    int KK = 0;
    NP++;

    for (int I = 0; I < NNQ; ++I) {
      if (I == NRESP - 1)
        continue;

      XV(KK, 0) = 1.0;
      XV(KK, 1) = A[I];
      SSE += std::pow(((A[I] - ZZSUM) / BTSUM) - ZTSUM[KK], 2);
      KK++;
    }

    // Compute YV matrix
    for (int K = 0; K < 2; ++K) {
      for (int J = 0; J < 2; ++J) {
        double SUM = 0.0;
        for (int I = 0; I < NQ; ++I) {
          SUM += XV(I, K) * XV(I, J);
        }
        YV(K, J) = SUM;
      }
    }

    // Compute inverse of YV
    double DET = 1.0 / (YV(0, 0) * YV(1, 1) - YV(0, 1) * YV(1, 0));
    double AA = YV(0, 0);
    double YV11 = YV(1, 1) * DET;
    double YV22 = AA * DET;
    double YV12 = -YV(0, 1) * DET;
    double YV21 = -YV(1, 0) * DET;

    YV(0, 0) = YV11;
    YV(1, 1) = YV22;
    YV(0, 1) = YV12;
    YV(1, 0) = YV21;

    // Compute ZV matrix
    for (int K = 0; K < NQ; ++K) {
      for (int I = 0; I < 2; ++I) {
        double SUM = 0.0;
        for (int J = 0; J < 2; ++J) {
          SUM += YV(I, J) * XV(K, J);
        }
        ZV(I, K) = SUM;
      }
    }

    // Compute A vector
    for (int I = 0; I < 2; ++I) {
      double SUM = 0.0;
      for (int J = 0; J < NQ; ++J) {
        SUM += ZV(I, J) * ZZ[J];
      }
      A[I] = SUM;
    }

    // Compute statistics
    double SUMA = 0.0;
    double SUMB = 0.0;
    double SUMC = 0.0;
    double SUMD = 0.0;
    double SUME = 0.0;
    for (int J = 0; J < NQ; ++J) {
      double AA = A[0] + A[1] * XV(J, 1);
      SUMA += AA * ZZ[J];
      SUMB += AA;
      SUMC += ZZ[J];
      SUMD += AA * AA;
      SUME += ZZ[J] * ZZ[J];
    }
    double aa = XNQ * SUMA - SUMB * SUMC;
    double bb = XNQ * SUMD - SUMB * SUMB;
    double cc = XNQ * SUME - SUMC * SUMC;
    double rsqrt = aa / std::sqrt(std::abs(bb * cc));
    if (A[1] < 0.0)
      rsqrt = -rsqrt;
    double RR = (aa * aa) / (bb * cc);

    double AA_result = A[0] + A[1] * XRES1;
    TSUM += AA_result;
    TTSUM += AA_result * AA_result;

    if (NRESP != 0)
      Q[NQ + NP - 1] = AA_result;

    // Transfer respondent coordinates
    int idx = (KDATA(IK - 1, 0) - 1) * 4;
    PSIMATRIX[idx] = A[0];
    PSIMATRIX[idx + 1] = A[1];
    PSIMATRIX[idx + 2] = (NRESP != 0) ? AA_result : 0.0;
    PSIMATRIX[idx + 3] = RR;

    if (A[1] >= 0.0) {
      IPOS++;
      QQ[IPOS + NQ - 1] = AA_result;
    }
    if (A[1] < 0.0) {
      INEG++;
      QQQ[INEG + NQ - 1] = AA_result;
    }

    for (int II = 0; II < NQ; ++II) {
      double AA = A[0] + XV(II, 1) * A[1];
      XSUM[II] += AA;
      XSUM[II + NQ] += AA * AA;
      double XCAND = ZZ[II];
      ASUM += AA;
      BSUM += XCAND;
      CSUM += AA * AA;
      DSUM += XCAND * XCAND;
      ESUM += AA * XCAND;
    }
  }

  double XNP = NP;
  SSE = std::sqrt(SSE / (XNP * XNQ));
  YSUM[NQ] = std::sqrt((TTSUM - (TSUM * TSUM) / XNP) / XNP);
  for (int I = 0; I < NQ; ++I) {
    ZSUM[I] = XSUM[I] / XNP;
    YSUM[I] = std::sqrt((XSUM[I + NQ] - (XSUM[I] * XSUM[I]) / XNP) / XNP);
  }
  ZSUM[NQ] = TSUM / XNP;

  // Transfer case numbers
  FITS[2] = NP;
  FITS[3] = IPOS;
  FITS[4] = INEG;

  // Compute R-square
  double RA = XNP * ESUM - ASUM * BSUM;
  double RB = XNP * CSUM - ASUM * ASUM;
  double RC = XNP * DSUM - BSUM * BSUM;
  double R1 = RA / std::sqrt(RB * RC);
  double R1S = R1 * R1;

  // Transfer R-square
  FITS[1] = R1S;

  // No need to deallocate in C++, as vectors and matrices are automatically managed
}
