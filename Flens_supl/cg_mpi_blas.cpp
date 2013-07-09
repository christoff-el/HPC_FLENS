#ifndef CG_MPI_BLAS_CPP
#define CG_MPI_BLAS_CPP 1

#include "cg_mpi_blas.h"

#include <iostream>

//FLENS-based MPI CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_mpi_blas(const MA &A, VX &x, const VB &b, VBC &bc,
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


template <typename MA, typename VX, typename VB, typename VBC>
int
cg_mpi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, VBC &fk_bc, 
						int maxIt, double tol)
{
    
    typedef int                                              IndexType;
    typedef flens::IndexOptions<IndexType, 1>         		 IndexBase;
    typedef flens::DenseVector<flens::Array<double> >		 DenseVector;
    
    //Check if sizes of matrices & vectors fit:
    assert(fk_A.numRows()==fk_b.values.length() && fk_A.numCols()==fk_x.values.length());
    assert(fk_x.type==typeI && fk_b.type==typeII);
    
    //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
	flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
	funk2flens_CRSmat(fk_A, fl_A);

	//Convert Funken DataVector b --> FLENS DenseVector:
	flens::FLENSDataVector fl_b(fk_b.values.length(), fk_b.coupling, (flens::VectorType)fk_b.type);
	funk2flens_DataVector(fk_b, fl_b);
		
	/***Solve using the FLENS-based CG solver ***/
	int iterCount;
	
	//x needs no MPI functionality, but we need attributes to copy to r1:
	flens::FLENSDataVector fl_x(fk_x.values.length(), fk_x.coupling, (flens::VectorType)fk_x.type);
	iterCount = cg_mpi_blas(fl_A, fl_x, fl_b, fk_bc, maxIt, tol);

	//Convert solution FLENSDataVector x --> Funken DataVector:
	flens2funk_DataVector(fl_x, fk_x);
	
	return iterCount;
}


#endif	//CG_MPI_BLAS_CPP