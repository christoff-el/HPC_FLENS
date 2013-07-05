#ifndef FLENS_DATA_VECTOR_H
#define FLENS_DATA_VECTOR_H 1

#include <flens/flens.cxx>

#include "flens2c.h"


namespace flens{

enum VectorType {typeI, typeII, nonMPI};

struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	//MPI Version (with Coupling, requires VectorType):
	explicit
    FLENSDataVector(int n, const Coupling &_coupling, const VectorType _vType)
        : 	DenseVector<Array<double> >(n),
        	vType(_vType),
        	coupling(_coupling)
    {
    }
    
    //Non-MPI version (no coupling, assumes 'nonMPI' VectorType):
    explicit
    FLENSDataVector(int n)
    	:	DenseVector<Array<double> >(n)
    {
    	vType = nonMPI;
    }

	//Member objects:
    VectorType vType;
	const Coupling &coupling;

	//Member methods:
	void typeII_2_I();
	void typeI_2_II();
	
	void commCrossPoints();
	void commBoundaryNodes();
	

}


void
FLENSDataVector::typeII_2_I()
{

	assert(type==typeII);
	
	//Sum up values at cross points:
	commCrossPoints();
	
	//Sum up values at boundary nodes:
	commBoundaryNodes();
	
	//Update VectorType:
	vType = typeI;

}

void
FLENSDataVector::typeI_2_II()
{

	assert(vType == typeI);
	
	//Divide values at cross points by the number of processes:
	for (int i=0; i<coupling.local2globalCrossPoints.length(); ++i) {
		(*this)(i+1) /= coupling.crossPointsNumProcs(i);
	}
	
	//Divide values at boundary nodes by 2 (since here 2 processes overlap):		(don't divide for cross points!)
	for (int i=0; i<coupling.numCoupling; ++i) {
		for (int j=1; j<coupling.boundaryNodes[i].length()-1; ++j) {
			(*this)(coupling.boundaryNodes[i](j)) /= 2.;
		}
	}
	
	//Update type:
	vType = typeII;

}

void
FLENSDataVector::commCrossPoints()
{
	
	DenseVector<Array<double> > u_crossPoints(coupling.numCrossPoints);
	
	//Local values at all global cross points:
	for (int i=0; i<coupling.local2globalCrossPoints.length(); ++i) {
		u_crossPoints(coupling.local2globalCrossPoints(i)) = (*this)(i+1);
	}
	
	double *u_crossPoints_tr;
	flens2c_DataVector(u_crossPoints, u_crossPoints_tr);
	
	double *u_crossPoints_gl = new double[coupling.numCrossPoints];
	 
	/*** MPI Communication for global cross points ***/
	MPI::COMM_WORLD.Allreduce(u_crossPoints_tr, u_crossPoints_gl, coupling.numCrossPoints,
										MPI::DOUBLE, MPI::SUM);
	
	for (int i=1; i<=coupling.local2globalCrossPoints.length(); ++i) {
		(*this)(i) = u_crossPoints_gl[coupling.local2globalCrossPoints(i)-1];
	}
	
	delete[] u_crossPoints_tr;
	delete[] u_crossPoints_gl;

}

void FLENSDataVector::commBoundaryNodes();
{

	for (int i=1; i<=coupling.maxColor; ++i) {
		for (int j=0; j<coupling.numCoupling; ++j) {
			if (coupling.colors(j) == i && coupling.boundaryNodes[j].length()-2 > 0) {
			
				//Only communicate if there is a boundary node on coupling boundary (no cross points):
				int sendLength = coupling.boundaryNodes[j].length()-2;
				
				DenseVector<Array<int> > sendIndex(sendLength);
				for (k=1; k<=sendLength; ++k) {					//copy manually until we can use blas::copy
					sendIndex(k) = coupling.boundaryNodes[j](k);
				}
				
				double *u_send = new double[sendLength];
				double *u_recv = new double[sendLength];
				
				//Set local values:
				for (int k=0; k<sendLength; ++k) {
					u_send[k] = (*this)(sendIndex(k));
				}
				
				//Get values from other processes:
				MPI::COMM_WORLD.Sendrecv(u_send, sendLength, MPI::DOUBLE, coupling.neighbourProcs(j)-1, 0,
											u_recv, sendLength, MPI::DOUBLE, coupling.neighbourProcs(j)-1,0);
				
				//Add values from other processes (!!numbering is opposite):
				for (int k=0; k<sendLength; ++k) {
					(*this)(sendIndex(k)) += u_recv[sendLength-k-1];
				}
				
			}
		}
	}

}





}	//namespace flens


#endif	//FLENS_DATA_VECTOR_H