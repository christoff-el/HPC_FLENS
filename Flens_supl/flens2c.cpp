/*#ifndef FLENS2C_CPP
#define FLENS2C_CPP 1

#include "flens2c.h"


//FLENS DenseVector --> C array:
void
flens2c_DataVector(flens::FLENSDataVector &fl_x, double *fk_x)
{

	fk_x = new double[fl_x.length()];
	
	//Copy FLENS vector to C array:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x[i-1] = fl_x(i);
	
	}

}


#endif	//FLENS2FUNK_CPP