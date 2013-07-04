#ifndef FLENS2FUNK_CPP
#define FLENS2FUNK_CPP 1

#include "flens2funk.h"


//FLENS DenseVector --> Funken Vector:
template <typename FLENSX, typename FUNKX>
void
flens2funk_Vector(FLENSX &fl_x, FUNKX &fk_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy FLENS vector to Funken vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x(i-1) = fl_x(i);
	
	}

}


#endif	//FLENS2FUNK_CPP