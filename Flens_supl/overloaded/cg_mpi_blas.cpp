#ifndef CG_MPI_BLAS_CPP
#define CG_MPI_BLAS_CPP 1

#include "cg_mpi_blas.h"


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
    typedef VX							VectorTypeI;
    typedef VB				    		VectorTypeII; 

	//Initialise variables:
    ElementType  	alpha, beta, rdot, rdotOld;
    VectorTypeI   	Ap(x), r1(x);
    VectorTypeII	p(b),  r2(b);	

    const ElementType  Zero(0), One(1);
    
    //r1 and r2 are the typeI and typeII vectors for the residuals


    //Set x to 0 at dirichlet nodes (where the value is fixed):
    for(int i=1; i<=bc.length(); ++i) {
          x(bc(i)) = Zero;
    }

    //Compute the residuals:
    //  p = A*x               ->  p is typeII
    //  r2 = b - A*x = b - p  -> r2 is typeII
    
    //blas::mv(NoTrans, One, A, x, Zero, p);
    p = A*x;
    
    //blas::axpy(-One, p, r2);
    r2 -= p;
    
    //Set r2 to 0 at the dirichlet nodes:
    for(int i=1; i<=bc.length(); ++i) {
    	r2(bc(i)) = Zero;
    }
    
    /*** MPI: Initialise direction p (as typeI residual ) ***/
    
    //Copy r2 to r1 - invokes II->I conversion:
    //blas::copy(r2, r1);
    r1 = r2;
    
    //Copy r1 to p - I -> I uses standard copy:
    //blas::copy(r1, p);
    p = r1;
    
    //Compute squared Norm of residuals, rdot = r*r:
    //rdot = blas::dot(r1, r2);
	rdot = r1*r2;
	
    for (int k=0; k<maxIt; k++) {
    
    	/*** Abort criterion ***/
        if (sqrt(rdot) <= tol) {
            return k;
        }

        //Compute  Ap = A*p, and set Ap to zero at dirichlet nodes:
        //blas::mv(NoTrans, One, A, p, Zero, Ap);
        Ap = A*p;

        for(int i=1; i<=bc.length(); ++i) {
        	Ap(bc(i)) = Zero;
        }
        
        //Compute alpha = rdot/(p * Ap):
        alpha = rdot; 
        
        //alpha /= blas::dot(Ap, p);
        alpha /= (Ap*p);
        
        //Update solution x by x += alpha*p:
        //blas::axpy(alpha, p, x);
        x += alpha*p;
        
        //Update local (= typeII) residual by r2 -= alpha*Ap:
        //blas::axpy(-alpha, Ap, r2);
        r2 -= alpha*Ap;

        /*** MPI: Get global (=type I) residual ***/
        //Copy r2 to r1 - invokes II->I conversion:
        //blas::copy(r2, r1);
        r1 = r2;
        
        //Compute  squared Norm of updated residuum rdot = r*r:
        rdotOld = rdot;
        
        //rdot = blas::dot(r1, r2);
        rdot = r1*r2;
        
        beta = rdot/rdotOld;
        
        //Update (global) direction p by p = beta*p + r1:
        //blas::scal(beta, p);
        //blas::axpy(One, r1, p);
        p = beta*p + r1;
        
    }
    
	return maxIt;    
}


#endif	//CG_MPI_BLAS_CPP