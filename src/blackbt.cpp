#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::depends(RcppEigen)]]

// The BLACKBT function using RcppEigen
// [[Rcpp::export]]
void BLACKBT(
    int NRESPONDENTS,
    int NISSUES,
    int NDIMENSIONS,
    int NP,
    int NY,
    int NF,
    int NFX,
    const Eigen::MatrixXd& XBIGONE,  // XBIGONE(NISSUES, NRESPONDENTS)
    Eigen::MatrixXd& XDATA,          // XDATA(NISSUES, NDIMENSIONS)
    Eigen::MatrixXd& W,              // W(NRESPONDENTS, NDIMENSIONS+2)
    double& SVSUM,
    Eigen::VectorXd& FITS2           // FITS2(6)
) {
  // Initialize variables
  int NF1 = NF + 1;
  Eigen::VectorXd DC(NY);   // Column means
  Eigen::VectorXd CC(NY);   // Column sum of squares
  std::vector<int> LL(NY, 0); // Non-missing counts for each column
  int KTOT = 0;
  SVSUM = 0.0;
  double SWSUM = 0.0;

  Eigen::MatrixXd X(NP, NY);
  Eigen::MatrixXd XS(NP, NY);
  Eigen::MatrixXd XSS(NP, NY);
  Eigen::MatrixXd XT(NP, NY);
  Eigen::MatrixXd PSIX(NP, NDIMENSIONS + 1);
  Eigen::MatrixXd XSAVE(NP, NDIMENSIONS + 2);
  Eigen::MatrixXd WSAVE(NY, NDIMENSIONS + 1);

  // Initialize DC and CC
  DC.setZero();
  CC.setZero();

  // Process data to compute sums and counts
  for (int I = 0; I < NP; ++I) {
    double SUM = 0.0;
    double SUMA = 0.0;
    int KK = 0;
    for (int J = 0; J < NY; ++J) {
      double XIJ = XBIGONE(I, J);
      X(I, J) = XIJ;
      if (std::abs(XIJ + 999.0) > 0.001) {
        SUM += XIJ * XIJ;
        SUMA += XIJ;
        KK++;
      }
      XS(I, J) = XIJ;
    }
    KTOT += KK;
    SVSUM += SUM;
    SWSUM += SUMA;
    for (int J = 0; J < NY; ++J) {
      double XIJ = X(I, J);
      XSS(I, J) = XIJ;
      if (std::abs(XIJ + 999.0) > 0.001) {
        LL[J]++;
        DC[J] += XIJ;
        CC[J] += XIJ * XIJ;
      }
    }
  }

  // Compute column means and sum of squares
  SVSUM = SVSUM - (SWSUM * SWSUM) / static_cast<double>(KTOT);
  int LTOT = NP * NY;
  double PMISS = static_cast<double>(LTOT - KTOT) / static_cast<double>(LTOT);
  LTOT = LTOT - KTOT;

  if (NF == 1) {
    FITS2(0) = NP;
    FITS2(1) = NY;
    FITS2(2) = KTOT;
    FITS2(3) = LTOT;
    FITS2(4) = PMISS * 100.0;
    FITS2(5) = SVSUM;
  }

  for (int J = 0; J < NY; ++J) {
    DC[J] = DC[J] / static_cast<double>(LL[J]);
    CC[J] = CC[J] - (DC[J] * DC[J]) * static_cast<double>(LL[J]);
  }

  // Master loop for NF dimensions
  for (int JJJ = 0; JJJ < NF; ++JJJ) {
    int KKK = JJJ;
    Eigen::VectorXd XX(2 * NY);
    XX.head(NY).setOnes();
    XX.segment(NY, NY).setZero();

    if (JJJ > 0) {
      XX.segment(NY, NY).setZero();
    } else {
      for (int J = 0; J < NY; ++J) {
        XX(NY + J) = -DC[J];
      }
    }

    double XXK = 0.0;
    double TXB = 0.0;
    int KKT = 0;

    // Compute initial estimates of PSIX
    for (int I = 0; I < NP; ++I) {
      int KK = 0;
      double SUM = 0.0;
      for (int J = 0; J < NY; ++J) {
        double XIJ = X(I, J);
        if (std::abs(XIJ + 999.0) <= 0.001) continue;
        SUM += (XIJ * XX[J] + XX[NY + J]) * static_cast<double>(LL[J]);
        KK++;
      }
      XT(I, 0) = -999.0;
      if (KK < NFX) continue;
      double WXB = SUM / static_cast<double>(KK);
      KKT++;
      TXB += WXB;
      XT(I, 0) = WXB;
      XXK += WXB * WXB;
    }

    TXB = TXB / static_cast<double>(KKT);
    XXK = XXK - static_cast<double>(KKT) * TXB * TXB;

    for (int I = 0; I < NP; ++I) {
      if (std::abs(XT(I, 0) + 999.0) <= 0.001) continue;
      XT(I, 0) = XT(I, 0) - TXB;
      XT(I, 1) = 1.0;
    }

    // Alternating Least Squares (ALS) loop
    double AREG = 0.0, BREG = 0.0;
    for (int MM = 0; MM < 4; ++MM) {
      // Call REGT (to be implemented)
      REGT(NRESPONDENTS, NISSUES, NDIMENSIONS, NP, 1, NY, W, XS, X, XT, KKK, AREG);

      // Call REG2T (to be implemented)
      double PXB = 0.0, PXS = 0.0;
      REG2T(NRESPONDENTS, NISSUES, NDIMENSIONS, NP, 1, NY, W, XS, X, XT, PXB, PXS, KKK, BREG);

      double QXS = PXS;
      double XCOR = std::sqrt(XXK / PXS);

      for (int I = 0; I < NP; ++I) {
        if (std::abs(XT(I, 0) + 999.0) <= 0.001) {
          PSIX(I, JJJ) = -999.0;
          continue;
        }
        PSIX(I, JJJ) = (XT(I, 0) - PXB) * XCOR;
        XT(I, 0) = (XT(I, 0) - PXB) * XCOR;
        double AA = PSIX(I, JJJ);
        if (AA <= -99.0) {
          XT(I, 0) = -999.0;
          PSIX(I, JJJ) = -999.0;
        }
      }
    }

    // Update X and XS matrices
    for (int I = 0; I < NP; ++I) {
      if (KKK == NF) {
        PSIX(I, NF) = 1.0;
      }
      for (int J = 0; J < NY; ++J) {
        if (KKK == NF)
          X(I, J) = XSS(I, J);
        XS(I, J) = X(I, J);
      }
    }

    if (KKK == NF)
      break;

    // Call CORR3 (to be implemented)
    CORR3(NRESPONDENTS, NISSUES, NP, NY, X, LL);
  }

  // Alternating Least Squares with full P matrix until convergence
  double AREG = 0.0, BREG = 0.0;
  for (int NN = 0; NN < 5; ++NN) {
    // Call REGT (to be implemented)
    REGT(NRESPONDENTS, NISSUES, NDIMENSIONS, NP, NF, NY, W, XS, X, PSIX, NF, AREG);

    // Call REG2T (to be implemented)
    double PXB = 0.0, PXS = 0.0;
    REG2T(NRESPONDENTS, NISSUES, NDIMENSIONS, NP, NF, NY, W, XS, X, PSIX, PXB, PXS, NF, BREG);

    // Center the PSI matrix at zero
    for (int K = 0; K < NF; ++K) {
      double SUM = PSIX.col(K).sum();
      SUM = SUM / static_cast<double>(NP);
      PSIX.col(K).array() -= SUM;
    }

    double AA = std::abs(AREG - BREG);
    if (AA < 0.01)
      break;
  }

  int NFPS = NF + 1;

  // Compute PW' + Jc' and put into X(NP, NY)
  for (int I = 0; I < NP; ++I) {
    for (int K = 0; K < NY; ++K) {
      XT(I, K) = 0.0;
    }
    if (std::abs(PSIX(I, 0) + 999.0) <= 0.001)
      continue;
    for (int K = 0; K < NY; ++K) {
      double SUM = 0.0;
      for (int J = 0; J < NFPS; ++J) {
        SUM += PSIX(I, J) * W(K, J);
      }
      X(I, K) = SUM;
    }
    for (int JJ = 0; JJ < NY; ++JJ) {
      if (std::abs(XS(I, JJ) + 999.0) > 0.001)
        XT(I, JJ) = XS(I, JJ);
      else
        XT(I, JJ) = X(I, JJ);
    }
  }

  // Singular Value Decomposition (SVD)
  Eigen::JacobiSVD<MatrixXd> svd1(XT, Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXd D = svd1.singularValues();

  // Compute residuals and sums
  double ESUM = 0.0;
  double SUMM = 0.0;
  double SUMM1 = 0.0;
  double SUMM2 = 0.0;

  for (int I = 0; I < NP; ++I) {
    for (int J = 0; J < NY; ++J) {
      SUMM += X(I, J) * X(I, J);
      if (std::abs(XS(I, J) + 999.0) > 0.001) {
        ESUM += std::pow(X(I, J) - XS(I, J), 2);
        SUMM1 += X(I, J) * X(I, J);
      } else {
        SUMM2 += X(I, J) * X(I, J);
      }
      // Store PW' + Jc' in XSS and PW' in XT
      XSS(I, J) = X(I, J);
      XT(I, J) = X(I, J) - W(J, NF);
    }
  }

  // SVD of PW' + Jc'
  Eigen::JacobiSVD<MatrixXd> svd2(XSS, Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXd DDD = svd2.singularValues();

  // SVD of PW'
  Eigen::JacobiSVD<MatrixXd> svd3(XT, Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXd DD = svd3.singularValues();

  // Update PSIX and W matrices
  MatrixXd U = svd3.matrixU();
  MatrixXd V = svd3.matrixV();
  VectorXd S = svd3.singularValues();

  for (int I = 0; I < NP; ++I) {
    for (int JJ = 0; JJ < NF; ++JJ) {
      PSIX(I, JJ) = V(I, JJ) * std::sqrt(S(JJ));
    }
  }

  for (int I = 0; I < NY; ++I) {
    for (int JJ = 0; JJ < NF; ++JJ) {
      W(I, JJ) = U(I, JJ) * std::sqrt(S(JJ));
    }
  }

  // Save PSI and W matrices
  XSAVE.leftCols(NF) = PSIX.leftCols(NF);
  WSAVE.leftCols(NF + 1) = W.leftCols(NF + 1);

  // Restore PSI and W matrices
  XDATA.leftCols(NF) = XSAVE.leftCols(NF);
  W.leftCols(NF + 1) = WSAVE.leftCols(NF + 1);
}
