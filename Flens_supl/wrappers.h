#ifndef WRAPPERS_H
#define WRAPPERS_H 1


#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FlensHeader.h"




//Wrapper: Funken --> FLENS --> Funken
int
cg_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol);

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
cg_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt, double tol);




//Wrapper: Funken --> FLENS --> Funken
int
gs_dense_nompi_blas_wrapper(Matrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							             int maxIt, double tol);

//Wrapper: Funken --> FLENS --> Funken
int
gs_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol);


//Wrapper: Funken --> FLENS >> GS >> FLENS --> Funken
int
gs_dense_mpi_blas_wrapper(Matrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt);

//Wrapper: Funken --> FLENS >> GS >> FLENS --> Funken
int
gs_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc,
						int maxIt);


#endif // WRAPPERS_H