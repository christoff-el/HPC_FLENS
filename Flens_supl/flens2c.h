/*#ifndef FLENS2C_H
#define FLENS2C_H 1

#include <flens/flens.cxx>

#include "FLENSDataVector.h"


//FLENS DataVector --> C array:
void
flens2c_DataVector(flens::FLENSDataVector &fl_x, double *fk_x);


#endif	//FLENS2C_H