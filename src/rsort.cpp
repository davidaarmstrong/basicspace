#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <cmath>

// Function to sort an array A and its associated index array IR
void RSORT(std::vector<double>& A, std::vector<int>& IR)
{
  int LA = A.size();
  if (LA <= 0) return;

  // Stack arrays for indices
  std::vector<int> IL(21);
  std::vector<int> IU(21);

  int M = 1;
  int I = 0;          // Adjusted for zero-based indexing
  int J = LA - 1;     // Adjusted for zero-based indexing
  double R = 0.375;

  begin_loop:
    if (I == J) goto next_segment;
    if (R > 0.5898437) goto adjust_R;
    R = R + 0.0390625;
    goto select_pivot;

    adjust_R:
      R = R - 0.21875;

    select_pivot:
      int K = I;
    // Select a central element of the array as pivot
    int IJ = I + static_cast<int>((J - I) * R);
    double T = A[IJ];
    int IT = IR[IJ];

    // First element of array is greater than T, interchange with T
    if (A[I] > T)
    {
      A[IJ] = A[I];
      A[I] = T;
      T = A[IJ];
      IR[IJ] = IR[I];
      IR[I] = IT;
      IT = IR[IJ];
    }

    int L = J;
    // If last element of array is less than T, interchange with T
    if (A[J] < T)
    {
      A[IJ] = A[J];
      A[J] = T;
      T = A[IJ];
      IR[IJ] = IR[J];
      IR[J] = IT;
      IT = IR[IJ];

      // If first element is greater than T, interchange with T
      if (A[I] > T)
      {
        A[IJ] = A[I];
        A[I] = T;
        T = A[IJ];
        IR[IJ] = IR[I];
        IR[I] = IT;
        IT = IR[IJ];
      }
    }

    partition_loop:
      // Find an element in the second half of the array which is smaller than T
      do {
      L = L - 1;
    } while (A[L] > T);

    // Find an element in the first half of the array which is greater than T
    do {
      K = K + 1;
    } while (A[K] < T);

    // Interchange these elements
    if (K <= L)
    {
      if (A[L] != A[K])
      {
        std::swap(A[L], A[K]);
        std::swap(IR[L], IR[K]);
      }
      goto partition_loop;
    }

    // Save upper and lower subscripts of the array yet to be sorted
    if ((L - I) <= (J - K))
    {
      IL[M - 1] = I;
      IU[M - 1] = L;
      I = K;
    }
    else
    {
      IL[M - 1] = K;
      IU[M - 1] = J;
      J = L;
    }
    M = M + 1;
    goto begin_loop;

    next_segment:
      M = M - 1;
    if (M == 0) return;
    I = IL[M - 1];
    J = IU[M - 1];

    if ((J - I) >= 11)
    {
      goto select_pivot;
    }

    if (I == 0) goto insertion_sort;

    I = I - 1;

    insertion_sort:
      while (true)
      {
        I = I + 1;
        if (I == J) goto next_segment;
        double T = A[I + 1];
        int IT = IR[I + 1];
        if (A[I] <= T) continue;
        K = I;
        do {
          A[K + 1] = A[K];
          IR[K + 1] = IR[K];
          K = K - 1;
        } while (K >= 0 && T < A[K]);
        A[K + 1] = T;
        IR[K + 1] = IT;
      }
}
