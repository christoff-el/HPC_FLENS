#ifndef FUNK2FLENS_H
#define FUNK2FLENS_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"

//Funken CRSMatrix --> FLENS CRS Matrix:
void
funk2flens_CRSmat(CRSMatrix &fk_A, flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1>> > &fl_A);

//Funken Vector --> FLENS DenseVector:
void
funk2flens_Vector(Vector &fk_x, flens::DenseVector<flens::Array<double> > &fl_x);


#endif	//FUNK2FLENS_H