#ifndef FLENS2FUNK_H
#define FLENS2FUNK_H 1

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"


//FLENS DenseVector --> Funken Vector:
template <typename FLENSX, typename FUNKX>
void
flens2funk_Vector(FLENSX &fl_x, FUNKX &fk_x);


#endif	//FLENS2FUNK_H