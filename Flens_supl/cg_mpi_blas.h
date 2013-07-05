#ifndef CG_MPI_BLAS_H
#define CG_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"

#include "funk2flens.h"
#include "flens2funk.h"


//FLENS-based CG solver:
template <typename MA, typename VX, typename VB>
int
cg_mpi_blas(const MA &A, const VB &b, VX &x,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon());
   

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
cg_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt, double tol);

#endif	//CG_MPI_BLAS_H