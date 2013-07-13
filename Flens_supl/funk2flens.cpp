#ifndef FUNK2FLENS_CPP
#define FUNK2FLENS_CPP 1

#include "funk2flens.h"


//Funken CRSMatrix --> FLENS CRS Matrix:
template <typename X=int>
void
funk2flens_CRSmat(CRSMatrix &fk_A, flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1>> > &fl_A)
{
	int n = fk_A.numRows();
	int m = fk_A.numCols();

	typedef int                                              				IndexType;
    typedef flens::IndexOptions<IndexType, 1>                				IndexBase;
    typedef flens::CoordStorage<double, flens::CoordRowColCmp, IndexBase>  	Coord;
	
	//FLENS CRS must be built from a Coordinate storage matrix:
    flens::GeCoordMatrix<Coord>  fl_A_coord(n, m);

	double *fk_data   = fk_A.data();
	int    *fk_colIdx = fk_A.colIndex();
	int    *fk_rowPtr = fk_A.rowPtr();
	
	
	for (int i=0; i<n; ++i) {
		for (int j=fk_rowPtr[i]; j<fk_rowPtr[i+1]; ++j) {
		
			fl_A_coord(i+1, fk_colIdx[j]+1) += fk_data[j];
			
		}
	}

	//Convert coordinate storage matrix to CRS:
	fl_A = fl_A_coord;

}

//Funken Matrix --> FLENS GeMatrix:
template <typename X=int>
void
funk2flens_mat(Matrix &fk_A, flens::GeMatrix<flens::FullStorage<double> > &fl_A)
{
	int n_lhs = fk_A.numRows(), m_lhs = fk_A.numCols();
	int n_rhs = fl_A.numRows(), m_rhs = fl_A.numCols();

	assert(n_lhs==n_rhs && m_lhs==m_rhs);

	// Copy Funken to FLENS
	for (int i=1; i<=n_rhs; ++i)
	{
		for (int j=1; j<=m_rhs; ++j)
		{
			fl_A(i,j) = fk_A(i-1,j-1);
		}
	}
};


//Funken Vector --> FLENS DenseVector:
template <typename X=int>
void
funk2flens_Vector(Vector &fk_x, flens::DenseVector<flens::Array<double> > &fl_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy Funken vector to FLENS vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fl_x(i) = fk_x(i-1);
	
	}

}

//Funken DataVector --> FLENSDataVector:
template <typename FLV>
void
funk2flens_DataVector(DataVector &fk_x, FLV &fl_x)
{
	
	assert(fl_x.length()==fk_x.values.length());
	
	//Add VectorType to FLENSDataVector (we must assume coupling already set, since it's const):
	//fl_x.vType = (flens::VectorType)fk_x.type;
	
	//Copy Funken data to FLENS data:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fl_x(i) = fk_x.values(i-1);
	
	}

}



#endif	//FUNK2FLENS_CPP