#ifndef FUNK2FLENS_H
#define FUNK2FLENS_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"

//Funken CRSMatrix --> FLENS CRS Matrix:
template <typename FUNKA, typename FLENSA>
void
funk2flens_CRSmat(FUNKA &fk_A, FLENSA &fl_A);

//Funken Vector --> FLENS DenseVector:
template <typename FUNKX, typename FLENSX>
void
funk2flens_Vector(FUNKX &fk_x, FLENSX &fl_x);


#endif	//FUNK2FLENS_H