#ifndef CG_MPI_BLAS_H
#define CG_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>


//FLENS-based CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon());


#include "cg_mpi_blas.cpp"

#endif	//CG_MPI_BLAS_H