#ifndef CG_MPI_BLAS_H
#define CG_MPI_BLAS_H 1

#include <flens/flens.cxx>
#include <limits>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"

#include "funk2flens.h"
#include "flens2funk.h"


//FLENS-based MPI CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   			int    maxIt = std::numeric_limits<int>::max(),
   			double tol = std::numeric_limits<double>::epsilon())
{
   	using namespace flens;

    typedef typename VB::ElementType  	ElementType;
    typedef typename VB::IndexType    	IndexType;
    typedef VB				    		VectorType;

	//Initialise variables:
    ElementType  	alpha, beta, rdot, rdotOld;
    VectorType   	Ap(x), p(x), r1(x), r2(b);	
    
    const ElementType  Zero(0), One(1);
    
    //r1 and r2 are the typeI and typeII vectors for the residuals

    
    //std::cout<< "b:"<<Ap.vType<<" x:"<<x.vType
    //Set x to 0 at dirichlet nodes (where the value is fixed):
    for(int i=0; i<bc.length(); ++i) {
          x(bc(i)) = 0;
    }
    
    //Compute the residuals:
    //  p = A*x               ->  p is typeII
    //  r2 = b - A*x = b - p  -> r2 is typeII
    blas::mv(NoTrans, One, A, x, Zero, p);
    blas::axpy(-One, p, r2);
    
    //Set r2 to 0 at the dirichlet nodes:
    for(int i=0; i<bc.length(); ++i) {
    	r2(bc(i)) = 0;
    }
    
    /*** MPI: Initialise direction p (as typeI residual ) ***/
    blas::copy(r2, r1);
    r1.typeII_2_I();
    //std::cout<<r1<<std::endl;
    blas::copy(r1, p);
    
    //Compute squared Norm of residuals, rdot = r*r:
    rdot = blas::dot(r1, r2);
    std::cout<<rdot<<std::endl;
    for (int k=0; k<maxIt; k++) {
    
    	/*** Abort criterion ***/
        if (sqrt(rdot) <= tol) {
            return k;
        }
        
        //Compute  Ap = A*p, and set Ap to zero at dirichlet nodes:
        blas::mv(NoTrans, One, A, p, Zero, Ap);

        for(int i=0; i<bc.length(); ++i) {
        	Ap(bc(i)) = 0;
        }
        
        //Compute alpha = rdot/(p * Ap):
        alpha = rdot; 
        alpha /= blas::dot(p, Ap);
        
        //Update solution x by x += alpha*p:
        blas::axpy(alpha, p, x);
        
        //Update (local = typeII) residual by r -= alpha*Ap:
        blas::axpy(-alpha, Ap, r2);

        /*** MPI: Get global (=type I) residual ***/
        blas::copy(r2, r1);
        r1.typeII_2_I();
        
        //Compute  squared Norm of updated residuum rdot = r*r:
        rdotOld = rdot;
        rdot = blas::dot(r1, r2);
        beta = rdot/rdotOld;
        
        //Update (global) direction p by p = beta*p + r1:
        blas::scal(beta, p);
        blas::axpy(One, r1, p);
        
    }
    
	return maxIt;    
}

	//Instatiations:
	template <>
	int cg_mpi_blas< flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1> > >,
								flens::FLENSDataVector, 
								flens::FLENSDataVector, 
								flens::DenseVector<flens::Array<int, flens::IndexOptions<int, 0> > > >
									(const flens::CRS<double, flens::IndexOptions<int, 1> > > &A,
										const flens::FLENSDataVector &b,
										flens::FLENSDataVector &x,
										flens::DenseVector<flens::Array<int, flens::IndexOptions<int, 0> > > > &bc,
										int maxIterations = std::numeric_limits<int>::max(),
										double tol = std::numeric_limits<double>::epsilon());
   

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
cg_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt, double tol);

#endif	//CG_MPI_BLAS_H