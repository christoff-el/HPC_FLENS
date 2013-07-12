#ifndef GS_MPI_BLAS_H
#define GS_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>
#include <mpi.h>
#include "../Fem/Coupling.hpp"
#include "FLENSDataVector.h"


//FLENS-based dense GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max());

//FLENS-based sparse GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max());
   

// TriDiag Solver
template <typename V>
void
solveTridiag(V &ldiag, V &diag, V &udiag, V &x, V &b);

#include  "gs_mpi_blas.cpp"

#endif	//GS_MPI_BLAS_H