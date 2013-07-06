#ifndef FLENS2FUNK_H
#define FLENS2FUNK_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "FLENSDataVector.h"

#include "../Fem/DataVector.hpp"


//FLENS DenseVector --> Funken Vector:
void
flens2funk_Vector(flens::DenseVector<flens::Array<double> > &fl_x, Vector &fk_x);

//FLENSDataVector --> Funken DataVector:
void
flens2funk_DataVector(flens::FLENSDataVector &fl_x, DataVector &fk_x);


#endif	//FLENS2FUNK_H