#ifndef CG_MPI_BLAS_CPP
#define CG_MPI_BLAS_CPP 1

#include "cg_mpi_blas.h"


//FLENS-based MPI CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon())
{
   	using namespace flens;

    typedef typename VB::ElementType  	ElementType;
    typedef typename VB::IndexType    	IndexType;
    typedef typename FLENSDataVector    VectorType;

	//Initialise variables:
    ElementType  	alpha, beta, rdot, rdotOld;
    VectorType   	Ap, p; r1, r2;	
    
    const ElementType  Zero(0), One(1);
    
    //r1 and r2 are the typeI and typeII vectors for the residuals
    
    //Copy b into r2:
    blas::copy(b, r2);
    
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
    blas::copy(r1, p);
    
    //Compute squared Norm of residuals, rdot = r*r:
    rdot = blas::dot(r, r);
    
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
        blas::dot(p, Ap, alpha);
        alpha = rdot / alpha;
        
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


int
cg_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						int maxIt, double tol)
{
    
    typedef int                                              IndexType;
    typedef flens::IndexOptions<IndexType, 1>         		 IndexBase;
    typedef flens::DenseVector<flens::Array<double> >		 DenseVector;
    
    //Check if sizes of matrices & vectors fit:
    assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
    assert(fk_x.type==typeI && fk_b.type==typeII);
    
    //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
	flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
	funk2flens_CRSmat(fk_A, fl_A);

	//Convert Funken DataVector b --> FLENS DenseVector:
	DenseVector fl_b(fk_b.length());
	funk2flens_Vector(fk_b, fl_b);
		
	//Solve using the FLENS-based CG solver:
	int iterCount;
	DenseVector fl_x(fl_b.length());
	iterCount = cg_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

	//Convert solution FLENS DenseVector x --> Funken Vector:
	flens2funk_Vector(fl_x, fk_x);
	
	return iterCount;
    

}


#endif	//CG_MPI_BLAS_CPP