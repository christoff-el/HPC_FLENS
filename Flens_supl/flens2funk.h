#ifndef FLENS2FUNK_H
#define FLENS2FUNK_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"


//FLENS DenseVector --> Funken Vector:
void
flens2funk_Vector(flens::DenseVector<flens::Array<double> > &fl_x, Vector &fk_x);


#endif	//FLENS2FUNK_H