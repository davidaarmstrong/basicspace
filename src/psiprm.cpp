#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

// Function to compute the Gram matrix R = PSI^T * PSI
void PSIPRM(
    int NP,
    int NF,
    const std::vector<std::vector<double>>& PSI,
    std::vector<std::vector<double>>& R,
    int IPRNT)
{
  // Suppress unused variable warning
  (void)IPRNT;

  // Initialize R matrix with zeros
  R.resize(NF, std::vector<double>(NF, 0.0));

  // Compute R(J,K) = sum over I of PSI[I][J] * PSI[I][K]
  for (int J = 0; J < NF; ++J) {
    for (int K = 0; K < NF; ++K) {
      double SUM = 0.0;
      for (int I = 0; I < NP; ++I) {
        SUM += PSI[I][J] * PSI[I][K];
      }
      R[J][K] = SUM;
    }
  }

  // Optional printing code is commented out in Fortran, so omitted here
}
