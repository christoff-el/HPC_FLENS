#ifndef FLENS2FUNK_H
#define FLENS2FUNK_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"

#include "../Fem/DataVector.hpp"


//FLENS DenseVector --> Funken IndexVector:
void
flens2funk_Vector(flens::DenseVector<flens::Array<int, flens::IndexOptions<int, 0> > > &fl_x, IndexVector &fk_x);

//FLENS DenseVector_Base0 --> Funken Vector:
void
flens2funk_Vector(flens::DenseVector<flens::Array<double, flens::IndexOptions<int, 0> > > &fl_x, Vector &fk_x);

//FLENS DenseVector --> Funken Vector:
void
flens2funk_Vector(flens::DenseVector<flens::Array<double> > &fl_x, Vector &fk_x);

//FLENSDataVector --> Funken DataVector:
void
flens2funk_DataVector(flens::FLENSDataVector &fl_x, DataVector &fk_x);


#endif	//FLENS2FUNK_H