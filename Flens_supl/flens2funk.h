#ifndef FLENS2FUNK_H
#define FLENS2FUNK_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"

#include "../Fem/DataVector.hpp"


//FLENS DenseVector --> Funken IndexVector:
template <typename FLV>
void
flens2funk_Vector(FLV &fl_x, IndexVector &fk_x);

//FLENS DenseVector_Base0 --> Funken Vector:
template <typename FLV>
void
flens2funk_Vector(flens::DenseVector<flens::Array<double, flens::IndexOptions<int, 0> > > &fl_x, Vector &fk_x);

//FLENS DenseVector --> Funken Vector:
template <typename FLV>
void
flens2funk_Vector(FLV &fl_x, Vector &fk_x);

//FLENSDataVector --> Funken DataVector:
template <typename FLV>
void
flens2funk_DataVector(FLV &fl_x, DataVector &fk_x);


#include "flens2funk.h"

#endif	//FLENS2FUNK_H