#include "Solver.hpp"

/********************************************************************************************************************/
/**********************************************     CG - Method      ************************************************/
/********************************************************************************************************************/

int CG(CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes,  int maxIt, double tol)
{
    
    //Solve for x using CG method under FLENS framework:
    maxIt = cg_nompi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);

    return maxIt;    
}

int CG_MPI(CRSMatrix &A, DataVector &x, DataVector &b, IndexVector &dirichletNodes,  int maxIt, double tol)
{
	
    //Solve for x using MPI CG method under FLENS framework:
    maxIt = cg_mpi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);
	
  	return maxIt;    
}

/********************************************************************************************************************/
/*******************************************   Gauß-Seidel - Method    **********************************************/
/********************************************************************************************************************/

int forwardGS( CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes, int maxIt, double tol)
{
    maxIt = gs_dense_nompi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);
        
    return maxIt;    
}

int forwardGS_MPI( CRSMatrix &A, DataVector &x, DataVector &b, Coupling &coupling, 
                   IndexVector &dirichletNodes, int maxIt)
{

    maxIt = gs_dense_mpi_blas_wrapper(A, x, b, coupling, dirichletNodes, maxIt);
        
    return maxIt;
}
