#ifndef CG_NOMPI_BLAS_H
#define CG_NOMPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"

#include "funk2flens.h"
#include "flens2funk.h"

#include "cg_nompi_blas.cpp"

//FLENS-based CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_nompi_blas(const MA &A, VX &x, const VB &b, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon());
   

//Wrapper: Funken --> FLENS --> Funken
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_nompi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, VBC &fk_bc,
							int maxIt, double tol);

#endif	//CG_NOMPI_BLAS_H