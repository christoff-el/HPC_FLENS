#ifndef CG_MPI_BLAS_CPP
#define CG_MPI_BLAS_CPP 1

#include "cg_mpi_blas.h"

#include <iostream>


int
cg_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
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
	iterCount = cg_mpi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

	//Convert solution FLENSDataVector x --> Funken DataVector:
	flens2funk_DataVector(fl_x, fk_x);
	
	return iterCount;
    

}


#endif	//CG_MPI_BLAS_CPP