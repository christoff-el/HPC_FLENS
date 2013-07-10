#ifndef GS_MPI_BLAS_H
#define GS_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"
#include "../Fem/Coupling.hpp"

#include "funk2flens.h"
#include "flens2funk.h"


//FLENS-based dense GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_mpi_blas(const MA &A, const VB &b, VX &x, Coupling &coupling, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max());

//FLENS-based sparse GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_mpi_blas(const MA &A, const VB &b, VX &x, Coupling &coupling, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max());
   

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
gs_dense_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt);

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
gs_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt);

template <typename V>
void
solveTridiag(V &ldiag, V &diag, V &udiag, V &x, V &b);

#endif	//GS_MPI_BLAS_H