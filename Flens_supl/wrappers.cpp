#ifndef WRAPPERS_CPP
#define WRAPPERS_CPP 1

#include "wrappers.h"


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
	iterCount = cg_nompi_blas(fl_A, fl_x, fl_b, fk_bc, maxIt, tol);

	//Convert solution FLENS DenseVector x --> Funken Vector:
	flens2funk_Vector(fl_x, fk_x);
	
	return iterCount;
};



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
	iterCount = cg_mpi_blas(fl_A, fl_x, fl_b, fk_bc, maxIt, tol);

	//Convert solution FLENSDataVector x --> Funken DataVector:
	flens2funk_DataVector(fl_x, fk_x);
	
	return iterCount;
};


//Wrapper: Funken --> FLENS --> Funken
int
gs_dense_nompi_blas_wrapper(Matrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							             int maxIt, double tol)
{
  typedef int                                          IndexType;
  typedef flens::IndexOptions<IndexType, 1>            IndexBase;
  typedef flens::DenseVector<flens::Array<double> >    DenseVector;
    
  //Check if sizes of matrices & vectors fit:
  assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
   
  //Convert Funken Matrix A --> FLENS GeMatrix:
  flens::GeMatrix<flens::FullStorage<double> > fl_A(fk_A.numRows(),fk_A.numCols());
  funk2flens_mat(fk_A, fl_A);

  //Convert Funken Vector b --> FLENS DenseVector:
  DenseVector fl_b(fk_b.length());
  funk2flens_Vector(fk_b, fl_b);
    
  //Solve using the FLENS-based GS solver:
  int iterCount;
  DenseVector fl_x(fl_b.length());
  iterCount = gs_dense_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

  //Convert solution FLENS DenseVector x --> Funken Vector:
  flens2funk_Vector(fl_x, fk_x);
  
  return iterCount;
};


//Wrapper: Funken --> FLENS --> Funken
int
gs_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
              int maxIt, double tol)
{
  typedef int                                          IndexType;
  typedef flens::IndexOptions<IndexType, 1>            IndexBase;
  typedef flens::DenseVector<flens::Array<double> >    DenseVector;
    
  //Check if sizes of matrices & vectors fit:
  assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
 
  //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
  flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
  funk2flens_CRSmat(fk_A, fl_A);

  //Convert Funken Vector b --> FLENS DenseVector:
  DenseVector fl_b(fk_b.length());
  funk2flens_Vector(fk_b, fl_b);
    
  //Solve using the FLENS-based GS solver:
  int iterCount;
  DenseVector fl_x(fl_b.length());
  iterCount = gs_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

  //Convert solution FLENS DenseVector x --> Funken Vector:
  flens2funk_Vector(fl_x, fk_x);
  
  return iterCount;
};


//Wrapper: Funken --> FLENS >> GS >> FLENS --> Funken
int
gs_dense_mpi_blas_wrapper(Matrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
						              int maxIt)
{
	typedef int                                              IndexType;
    typedef flens::IndexOptions<IndexType, 1>         		 IndexBase;
    typedef flens::DenseVector<flens::Array<double> >		 DenseVector;
    
    //Check if sizes of matrices & vectors fit:
    assert(fk_A.numRows()==fk_b.values.length() && fk_A.numCols()==fk_x.values.length());
    assert(fk_x.type==typeI && fk_b.type==typeII);
    
    //Convert Funken Matrix A --> FLENS GeMatrix:
	flens::GeMatrix<flens::FullStorage<double> > fl_A(fk_A.numRows(),fk_A.numCols());
	funk2flens_mat(fk_A, fl_A);

	//Convert Funken DataVector b --> FLENS DenseVector:
	flens::FLENSDataVector fl_b(fk_b.values.length(), fk_b.coupling, (flens::VectorType)fk_b.type);
	funk2flens_DataVector(fk_b, fl_b);
		
	/***Solve using the FLENS-based GS solver ***/
	int iterCount;
	
	//Convert Funken DataVector x --> FLENS DenseVector:
	flens::FLENSDataVector fl_x(fk_x.values.length(), fk_x.coupling, (flens::VectorType)fk_x.type);
	iterCount = gs_dense_mpi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt);

	//Convert solution FLENSDataVector x --> Funken DataVector:
	flens2funk_DataVector(fl_x, fk_x);
	
	return iterCount;
};

//Wrapper: Funken --> FLENS >> CG >> FLENS --> Funken
int
gs_mpi_blas_wrapper(CRSMatrix &fk_A, DataVector &fk_x, DataVector &fk_b, IndexVector &fk_bc, 
                        int maxIt)
{
    typedef int                                              IndexType;
    typedef flens::IndexOptions<IndexType, 1>                IndexBase;
    typedef flens::DenseVector<flens::Array<double> >        DenseVector;
    
    //Check if sizes of matrices & vectors fit:
    assert(fk_A.numRows()==fk_b.values.length() && fk_A.numCols()==fk_x.values.length());
    assert(fk_x.type==typeI && fk_b.type==typeII);
    
    //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
    flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
    funk2flens_CRSmat(fk_A, fl_A);

    //Convert Funken DataVector b --> FLENS DenseVector:
    flens::FLENSDataVector fl_b(fk_b.values.length(), fk_b.coupling, (flens::VectorType)fk_b.type);
    funk2flens_DataVector(fk_b, fl_b);
        
    /***Solve using the FLENS-based GS solver ***/
    int iterCount;
    
    //Convert Funken DataVector x --> FLENS DenseVector:
    flens::FLENSDataVector fl_x(fk_x.values.length(), fk_x.coupling, (flens::VectorType)fk_x.type);
    iterCount = gs_mpi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt);

    //Convert solution FLENSDataVector x --> Funken DataVector:
    flens2funk_DataVector(fl_x, fk_x);
    
    return iterCount;
};


#endif // WRAPPERS_CPP