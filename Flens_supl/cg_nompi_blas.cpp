#ifndef CG_NOMPI_BLAS_CPP
#define CG_NOMPI_BLAS_CPP 1

#include "cg_nompi_blas.h"




//Wrapper: Funken --> FLENS --> Funken
int
cg_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
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
	iterCount = cg_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

	//Convert solution FLENS DenseVector x --> Funken Vector:
	flens2funk_Vector(fl_x, fk_x);
	
	return iterCount;
}

#endif	//CG_NOMPI_BLAS_CPP