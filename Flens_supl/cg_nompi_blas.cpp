#ifndef CG_NOMPI_BLAS_CPP
#define CG_NOMPI_BLAS_CPP 1

#include "cg_nompi_blas.h"


//FLENS-based CG solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_nompi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon())
{
    using namespace flens;

    typedef typename VB::ElementType  ElementType;
    typedef typename VB::IndexType    IndexType;
    typedef typename VB::NoView       VectorType;

    ElementType  alpha, beta, rNormSquare, rNormSquarePrev;
    VectorType   Ap, r, p;

    const ElementType  Zero(0), One(1);

	//Set x to 0 at dirichlet nodes (where the value is fixed):
    for(int i=0; i<bc.length(); ++i) {
        x(bc(i)) = 0;
    }
    
    blas::copy(b, r);
    blas::mv(NoTrans, -One, A, x, One, r);

	//Set r to 0 at dirichlet nodes:
    for(int i=0; i<bc.length(); ++i) {
        r(bc(i)) = 0;
    }
    
    blas::copy(r, p);

    rNormSquare = blas::dot(r, r);

    for (int k=1; k<=maxIterations; ++k) {

        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }

        blas::mv(NoTrans, One, A, p, Zero, Ap);

 		//Set Ap to zero at dirichlet nodes:
        for(int i=0; i<bc.length(); ++i) {
        	Ap(bc(i)) = 0;
    	}
    	
        alpha = rNormSquare/blas::dot(p, Ap);

        blas::axpy(alpha, p, x);

        blas::axpy(-alpha, Ap, r);

        rNormSquarePrev = rNormSquare;

        rNormSquare = blas::dot(r, r);

        beta = rNormSquare/rNormSquarePrev;

        blas::scal(beta, p);
        blas::axpy(One, r, p);
        
    }
    
    return maxIterations;
}


//Wrapper: Funken --> FLENS --> Funken
template <typename MA, typename VX, typename VB, typename VBC>
int
cg_nompi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, VBC &fk_bc,
							int maxIt, double tol)
{

	typedef int                                              IndexType;
    typedef flens::IndexOptions<IndexType, 1>         		 IndexBase;
    typedef flens::DenseVector<flens::Array<double> >		 DenseVector;
    
    //Check if sizes of matrices & vectors fit:
    assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
    
	//Convert Funken CRSMatrix A --> FLENS CRS Matrix:
	flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
	funk2flens_CRSmat(fk_A, fl_A);

	//Convert Funken Vector b --> FLENS DenseVector:
	DenseVector fl_b(fk_b.length());
	funk2flens_Vector(fk_b, fl_b);
		
	//Solve using the FLENS-based CG solver:
	int iterCount;
	DenseVector fl_x(fl_b.length());
	iterCount = cg_nompi_blas(fl_A, fl_x, fl_b, fk_bc, maxIt, tol);

	//Convert solution FLENS DenseVector x --> Funken Vector:
	flens2funk_Vector(fl_x, fk_x);
	
	return iterCount;
}

#endif	//CG_NOMPI_BLAS_CPP