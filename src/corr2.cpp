#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <cmath>

// Function to perform the correlation calculations
// [[Rcpp::export]]
void CORR2(
    int NRESPONDENTS,
    int NISSUES,
    int NP,
    int NY,
    const std::vector<double>& X,
    std::vector<double>& R,
    std::vector<int>& LL,
    std::vector<int>& MPOS,
    int& KS,
    int& KPOS,
    int IPRNT)
{
  // Allocate and initialize SA, SB, SC, SD matrices
  std::vector<std::vector<double>> SA(NY, std::vector<double>(NY, 0.0));
  std::vector<std::vector<double>> SB(NY, std::vector<double>(NY, 0.0));
  std::vector<std::vector<double>> SC(NY, std::vector<double>(NY, 0.0));
  std::vector<std::vector<double>> SD(NY, std::vector<double>(NY, 0.0));

  // Initialize SA, SB, SC, SD to zero (already initialized during allocation)
  // Begin main computation
  for (int I = 0; I < NP; ++I) {
    for (int J = 0; J < NY; ++J) {
      // Skip if condition is met
      if (std::abs(X[I][J] + 999.0) <= 0.001) continue;

      bool skip_J = false;
      for (int JJ = 0; JJ <= J; ++JJ) {
        if (std::abs(X[I][JJ] + 999.0) <= 0.001) {
          skip_J = true;
          break;
        }

        // Update SA matrix
        SA[J][JJ] += X[I][J];
        if (J != JJ)
          SA[JJ][J] += X[I][JJ];

        // Update SB matrix
        SB[J][JJ] += X[I][J] * X[I][J];
        if (J != JJ)
          SB[JJ][J] += X[I][JJ] * X[I][JJ];

        // Update SC matrix
        SC[J][JJ] += X[I][J] * X[I][JJ];
        SC[JJ][J] = SC[J][JJ];

        // Update SD matrix
        SD[J][JJ] += 1.0;
      }
      if (skip_J) continue;
    }
  }

  // Compute correlation matrix R
  for (int J = 0; J < NY; ++J) {
    for (int JJ = 0; JJ <= J; ++JJ) {
      double AA = SD[J][JJ] * SC[J][JJ] - SA[J][JJ] * SA[JJ][J];
      double BB = SD[J][JJ] * SB[J][JJ] - SA[J][JJ] * SA[J][JJ];
      double CC = SD[J][JJ] * SB[JJ][J] - SA[JJ][J] * SA[JJ][J];
      double BBCC = BB * CC;

      if (BBCC <= 0.0) {
        R[JJ][J] = 0.0;
      } else {
        R[JJ][J] = AA / std::sqrt(BBCC);
      }
      R[J][JJ] = R[JJ][J];
    }
  }

  // Find the row with the largest total sum of absolute correlations
  double BB = -99.0;
  KS = 0;
  for (int J = 0; J < NY; ++J) {
    double SUM = 0.0;
    for (int JJ = 0; JJ < NY; ++JJ) {
      SUM += std::abs(R[J][JJ]);
    }
    if (SUM > BB) {
      BB = SUM;
      KS = J;
    }
  }
  int KSOLD = KS;

  // Initialize LL vector based on R matrix
  for (int J = 0; J < NY; ++J) {
    if (R[KS][J] <= 0.0)
      LL[J] = -1;
    else
      LL[J] = 1;
  }

  // Iteratively change signs of columns to maximize positive entries in R
  int NYD2 = (NY - 1) / 2;
  KPOS = 0;
  for (int JK = 0; JK < NY; ++JK) {
    for (int J = 0; J < NY; ++J) {
      int KK = 0;
      int KSUM = 0;
      for (int JJ = 0; JJ < NY; ++JJ) {
        double AA = R[J][JJ] * static_cast<double>(LL[JJ]) * static_cast<double>(LL[J]);
        if (JK == NY - 1 && AA >= 0.0) {
          KPOS++;
          KSUM++;
        }
        if (AA < 0.0)
          KK++;
      }
      if (KK > NYD2) {
        LL[J] = LL[J] * (-1);
        KS = 999;
      }
      if (JK == NY - 1)
        MPOS[J] = KSUM;
    }
  }

  // The SA, SB, SC, SD vectors are automatically deallocated when they go out of scope
}
\
