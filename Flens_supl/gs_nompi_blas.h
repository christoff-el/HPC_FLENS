#ifndef GS_NOMPI_BLAS_H
#define GS_NOMPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"

#include "funk2flens.h"
#include "flens2funk.h"

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
int
gs_dense_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol);

//Wrapper: Funken --> FLENS --> Funken
int
gs_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol);

#endif // GS_NOMPI_BLAS_H