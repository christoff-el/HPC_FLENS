#ifndef FLENS2FUNK_CPP
#define FLENS2FUNK_CPP 1

#include "flens2funk.h"


//FLENS DenseVector --> Funken IndexVector:
template <typename FLV>
void
flens2funk_Vector(FLV &fl_x, IndexVector &fk_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy FLENS vector to Funken vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x(i-1) = fl_x(i);
	
	}

}

//FLENS DenseVector_Base0 --> Funken Vector:
template <typename FLV>
void
flens2funk_Vector(flens::DenseVector<flens::Array<double, flens::IndexOptions<int, 0> > > &fl_x, Vector &fk_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy FLENS vector to Funken vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x(i-1) = fl_x(i);
	
	}

}

//FLENS DenseVector --> Funken Vector:
template <typename FLV>
void
flens2funk_Vector(FLV &fl_x, Vector &fk_x)
{

	assert(fl_x.length()==fk_x.length());
	
	//Copy FLENS vector to Funken vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x(i-1) = fl_x(i);
	
	}

}

//FLENSDataVector --> Funken DataVector:
template <typename FLV>
void
flens2funk_DataVector(FLV &fl_x, DataVector &fk_x)
{

	assert(fl_x.length()==fk_x.values.length());
	
	//Add VectorType to Funken DataVector (we must assume coupling already set, since it's const):
	//fk_x.type = (vectorType)fl_x.vType;

	//Copy FLENS vector to Funken vector:
	for (int i=1; i<=fl_x.length(); ++i) {
	
		fk_x.values(i-1) = fl_x(i);
	
	}

}


#endif	//FLENS2FUNK_CPP