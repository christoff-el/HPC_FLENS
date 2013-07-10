#ifndef GS_NOMPI_BLAS_H
#define GS_NOMPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FlensHeader.h"


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



//Wrapper: Funken --> FLENS --> Funken
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_nompi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, VBC &fk_bc,
							int maxIt, double tol);

<<<<<<< HEAD
//Wrapper: Funken --> FLENS --> Funken
int
gs_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol);
=======

#include "gs_nompi_blas.cpp"
>>>>>>> chris

#endif // GS_NOMPI_BLAS_H