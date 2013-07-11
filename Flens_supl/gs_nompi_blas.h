#ifndef GS_NOMPI_BLAS_H
#define GS_NOMPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>


//FLENS-based dense GS solver
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_nompi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon());

//FLENS-based sparse GS solver
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_nompi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon());

#include "gs_nompi_blas.cpp"

#endif // GS_NOMPI_BLAS_H