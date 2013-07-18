#ifndef CRS_SORT_H
#define CRS_SORT_H 1

#include <stdio.h>
#include <flens/flens.cxx>

typedef flens::DenseVector<flens::Array<double> > 	DVector;
typedef flens::DenseVector<flens::Array<int> > 		IVector;

namespace flens{

void
CRS_sort(IVector &_data, IVector &row, IVector &col, IVector &values);

void
CRS_insort(int* index, int* data, int length);

}; //namespace flens


#endif //CRS_SORT_H