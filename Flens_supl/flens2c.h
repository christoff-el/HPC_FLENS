#ifndef FLENS2C_H
#define FLENS2C_H 1

#include <flens/flens.cxx>


//FLENS DenseVector --> C array:
void
flens2c_DenseVector(flens::DenseVector<flens::Array<double> > &fl_x, double *fk_x);


#endif	//FLENS2C_H