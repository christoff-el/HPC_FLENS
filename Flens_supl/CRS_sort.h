#ifndef CRS_SORT_H
#define CRS_SORT_H 1

#include <stdio.h>
#include <flens/flens.cxx>

typedef flens::DenseVector<flens::Array<double> > 	DVector;
typedef flens::DenseVector<flens::Array<int> > 		IVector;

namespace flens{

template <typename DVector, typename IVector>
void
CRS_sort(DVector &_data, IVector &row, IVector &col, DVector &values);

template <typename Data, typename Index>
void
CRS_insort(Index* index, Data* data, Index length);

}; //namespace flens

#include "CRS_sort.cpp"

#endif //CRS_SORT_H
