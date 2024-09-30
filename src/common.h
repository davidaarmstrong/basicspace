#ifndef COMMON_H
#define COMMON_H

// Include necessary libraries
#include <RcppEigen.h>
#include <iostream>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

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
);

// Function prototypes
void BLACKB(
    int NP,
    int NRESPONDENTS,
    int NISSUES,
    int NDIMENSIONS,
    int NF,
    const Eigen::MatrixXd& XBIGONE,    // XBIGONE(NRESPONDENTS, NISSUES)
    Eigen::MatrixXd& XDATA,            // XDATA(NRESPONDENTS, NDIMENSIONS)
    Eigen::MatrixXd& W,                // W(NISSUES, NDIMENSIONS+2)
    double& SVSUM,
    Eigen::VectorXd& FITS2,            // FITS2(6)
    int IPRNT
);

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
);

void CORR2(int NRESPONDENTS,
           int NISSUES,
           int NP,
           int NY,
           Eigen::MatrixXd& X,         // X(NRESPONDENTS, NISSUES)
           Eigen::MatrixXd& R,         // R(NISSUES, NISSUES)
           Eigen::VectorXi& LL,        // LL(NISSUES)
           Eigen::VectorXi& MPOS,      // MPOS(NISSUES)
           int KS,
           int KPOS,
           int IPRNT);

void CORR3(
    int NRESPONDENTS,
    int NISSUES,
    int NP,
    int NY,
    const std::vector<std::vector<double>>& X,
    std::vector<int>& LL,
    std::vector<int>& MPOS,
    int& KS,
    int& KPOS);

void CORR22(
    int NRESPONDENTS,
    int NISSUES,
    int NP,
    int NY,
    const Eigen::MatrixXd& X,         // Input matrix X (Eigen::MatrixXd)
    Eigen::MatrixXd& R,               // Correlation matrix R (Eigen::MatrixXd)
    Eigen::VectorXi& LL,              // Integer vector LL
    Eigen::VectorXi& MPOS,            // Integer vector MPOS
    int& KS,
    int& KPOS,
    int IPRNT);

void PSIPRM(
    int NP,
    int NF,
    const std::vector<std::vector<double>>& PSI,
    std::vector<std::vector<double>>& R,
    int IPRNT);

void REG(int NRESPONDENTS,
         int NISSUES,
         int NDIMENSIONS,
         int NP,
         int NF,
         int NY,
         Eigen::VectorXd& TSUM,                // TSUM(2 * NISSUES)
         Eigen::MatrixXd& W,                   // W(NISSUES, NDIMENSIONS+2)
         Eigen::MatrixXd& XS,                  // XS(NRESPONDENTS, NISSUES)
         Eigen::MatrixXd& X,                   // X(NRESPONDENTS, NISSUES)
         Eigen::MatrixXd& PSI,                 // PSI(NRESPONDENTS, NISSUES)
         int IPRNT,
         int ILAST,
         int KKK,
         double AREG,
         double BREG);

void REG2(int NRESPONDENTS,
          int NISSUES,
          int NDIMENSIONS,
          int NP,
          int NF,
          int NY,
          Eigen::MatrixXd& W,  // W(NISSUES, NDIMENSIONS+2)
          Eigen::MatrixXd& XS, // XS(NRESPONDENTS, NISSUES)
          Eigen::MatrixXd& X,  // X(NRESPONDENTS, NISSUES)
          Eigen::MatrixXd& PSI,// PSI(NRESPONDENTS, NISSUES)
          double& PXB,
          double& PXS,
          int KKK,
          int IPRNT,
          double AREG,
          double& BREG);

void REG2T(int NRESPONDENTS,
           int NISSUES,
           int NDIMENSIONS,
           int NP,
           int NF,
           int NY,
           Eigen::MatrixXd& W,     // W(NISSUES, NDIMENSIONS+2)
           Eigen::MatrixXd& XS,    // XS(NISSUES, NRESPONDENTS)
           Eigen::MatrixXd& X,     // X(NISSUES, NRESPONDENTS)
           Eigen::MatrixXd& PSI,   // PSI(NISSUES, NRESPONDENTS)
           double& PXB,
           double& PXS,
           int KKK,
           int NWHO,
           double& BREG);

void REGA(int NISSUES,
          int NDIMENSIONS,
          int NS,
          int NF,
          const Eigen::MatrixXd& A,
          const Eigen::VectorXd& Y,
          Eigen::VectorXd& V,
          int IPRNT);

void REGAT(int NRESPONDENTS,
           int NDIMENSIONS,
           int NS,
           int NF,
           const Eigen::MatrixXd& A,
           const Eigen::VectorXd& Y,
           Eigen::VectorXd& V);

void REGT(int NRESPONDENTS,
          int NISSUES,
          int NDIMENSIONS,
          int NP,
          int NF,
          int NY,
          std::vector<double>& TSUM,    // TSUM(2*NRESPONDENTS)
          Eigen::MatrixXd& W,           // W(NRESPONDENTS, NDIMENSIONS+2)
          Eigen::MatrixXd& XS,          // XS(NISSUES, NRESPONDENTS)
          Eigen::MatrixXd& X,           // X(NISSUES, NRESPONDENTS)
          Eigen::MatrixXd& PSI,         // PSI(NISSUES, NRESPONDENTS)
          int IPRNT,
          int ILAST,
          int KKK,
          double& AREG);

void RSORT(std::vector<double>& A,
           std::vector<int>& IR);

void RSQUR(
    int NRESPONDENTS,
    int NISSUES,
    int NP,
    int NY,
    double& R,
    const Eigen::MatrixXd& A,  // Use Eigen::MatrixXd for A
    const Eigen::MatrixXd& B,  // Use Eigen::MatrixXd for B
    int IPRNT);




// Add other function declarations here...

#endif // COMMON_H

