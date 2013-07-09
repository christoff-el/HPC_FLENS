#ifndef GS_MPI_BLAS_H
#define GS_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>
#include <mpi.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FlensHeader.h"
#include "../Fem/Coupling.hpp"


//FLENS-based GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_mpi_blas(const MA &A, const VB &b, VX &x, Coupling &coupling, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max());
   

//Wrapper: Funken --> FLENS >> GS >> FLENS --> Funken
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_mpi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, Coupling &coupling, VBC &fk_bc, 
						int maxIt);

template <typename V>
void
solveTridiag(V &ldiag, V &diag, V &udiag, V &x, V &b);


#include "gs_mpi_blas.cpp"

#endif	//GS_MPI_BLAS_H