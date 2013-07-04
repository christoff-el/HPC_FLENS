#ifndef FUNK2FLENS_CPP
#define FUNK2FLENS_CPP 1

#include "funk2flens.h"


//Funken CRSMatrix --> FLENS CRS Matrix:
template <typename FUNKA, typename FLENSA>
void
funk2flens_CRSmat(const FUNKA &fk_A, FLENSA &fl_A)
{

	typedef int                                              				IndexType;
    typedef flens::IndexOptions<IndexType, 1>                				IndexBase;
    typedef flens::CoordStorage<double, flens::CoordRowColCmp, IndexBase>  	Coord;

	int n = fk_A.numRows();
	int m = fk_A.numCols();
	
    flens::GeCoordMatrix<Coord>  fl_A_coord(n, m);

	double *fk_data   = fk_A.data();
	int    *fk_colIdx = fk_A.colIndex();
	int    *fk_rowPtr = fk_A.rowPtr();
	
	for (int i=0; i<n; ++i) {
		for (int j=fk_rowPtr[i]; j<=fk_rowPtr[i+1]; ++j) {
		
			fl_A_coord(i+1, fk_colIdx[j]+1) = fk_data[j];
			
		}
	}
	
	fl_A = fl_A_coord;

}


//Funken Vector --> FLENS DenseVector:
template <typename FUNKX, typename FLENSX>
void
funk2flens_Vector(const FUNKX &fk_x, FLENSX &fl_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy Funken vector to FLENS vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fl_x(i) = fk_x(i-1);
	
	}

}



#endif	//FUNK2FLENS_CPP