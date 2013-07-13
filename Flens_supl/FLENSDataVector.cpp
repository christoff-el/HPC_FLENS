#ifndef FLENS_DATA_VECTOR_CPP
#define FLENS_DATA_VECTOR_CPP 1

#include "FLENSDataVector.h"

namespace flens{


//Non-MPI --> Restricted to only FLnonMPI specialisation.
template <typename VTYPE>
FLENSDataVector<VTYPE>::FLENSDataVector(int n)
{
	VTYPE::CHK;			//<-- If scope ever reaches here, compilation will fail.
						//		e.g. if FLENSDataVector<double> instantiated.
}

template <>
FLENSDataVector<FLvNonMPI>::FLENSDataVector(int n)
	: 	DenseVector<Array<double> >(n),
		coupling(NULL)
{
	//Permits instatiation of FlNonMPI specialisation.
}


//MPI --> Restricted to only FLvTypeI, FLvTypeII specialisations.
template <typename VTYPE>
FLENSDataVector<VTYPE>::FLENSDataVector(int n, const Coupling &_coupling)
{
	VTYPE::CHK;			//<-- If scope ever reaches here, compilation will fail.
						//		e.g. if FLENSDataVector<double> instantiated.
}

template <>
FLENSDataVector<FLvTypeI>::FLENSDataVector(int n, const Coupling &_coupling)
	:	DenseVector<Array<double> >(n),
		coupling(_coupling)
{
	//Permits instatiation of FlvTypeI specialisation.
}

template <>
FLENSDataVector<FLvTypeII>::FLENSDataVector(int n, const Coupling &_coupling)
	:	DenseVector<Array<double> >(n),
		coupling(_coupling)
{
	//Permits instatiation of FlvTypeII specialisation.
}
	
	

template <>	
void
FLENSDataVector<FLvTypeI>::typeII_2_I()
{

	//assert(vType==typeII);
	
	//Sum up values at cross points:
	commCrossPoints();

	//Sum up values at boundary nodes:
	commBoundaryNodes();
	
	//Update VectorType:
	//vType = typeI;

}

template <>
void
FLENSDataVector<FlvTypeI>::typeI_2_II()
{

	//assert(vType == typeI);
	
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

template <typename VTYPE>
void
FLENSDataVector<VTYPE>::commCrossPoints()
{
	
	FLENSDataVector u_crossPoints(coupling.numCrossPoints);

	//Local values at all global cross points:
	for (int i=0; i<coupling.local2globalCrossPoints.length(); ++i) {
		u_crossPoints(coupling.local2globalCrossPoints(i)) = (*this)(i+1);
	}

	double *u_crossPoints_tr = u_crossPoints.vec2c();
	
	//std::cout<<u_crossPoints_tr[1]<<std::endl;
	double *u_crossPoints_gl = new double[coupling.numCrossPoints];
	 
	/*** MPI Communication for global cross points ***/
	MPI::COMM_WORLD.Allreduce(u_crossPoints_tr, u_crossPoints_gl, coupling.numCrossPoints,
										MPI::DOUBLE, MPI::SUM);

	for (int i=0; i<coupling.local2globalCrossPoints.length(); ++i) {
		(*this)(i+1) = u_crossPoints_gl[coupling.local2globalCrossPoints(i)-1];
	}

	delete[] u_crossPoints_tr;
	delete[] u_crossPoints_gl;

}

template <typename VTYPE>
void 
FLENSDataVector<VTYPE>::commBoundaryNodes()
{

	for (int i=1; i<=coupling.maxColor; ++i) {
		for (int j=0; j<coupling.numCoupling; ++j) {
			if (coupling.colors(j) == i && coupling.boundaryNodes[j].length()-2 > 0) {
			
				//Only communicate if there is a boundary node on coupling boundary (no cross points):
				int sendLength = coupling.boundaryNodes[j].length()-2;

				DenseVector<Array<int> > sendIndex(sendLength);
				for (int k=0; k<sendLength; ++k) {					//copy manually until we can use blas::copy
					sendIndex(k+1) = coupling.boundaryNodes[j](k+1);
				}

				double *u_send = new double[sendLength];
				double *u_recv = new double[sendLength];
				
				//Set local values:
				for (int k=0; k<sendLength; ++k) {
					u_send[k] = (*this)(sendIndex(k+1));
				}

				//Get values from other processes:
				MPI::COMM_WORLD.Sendrecv(u_send, sendLength, MPI::DOUBLE, coupling.neighbourProcs(j)-1, 0,
											u_recv, sendLength, MPI::DOUBLE, coupling.neighbourProcs(j)-1,0);
				
				//Add values from other processes (!!numbering is opposite):
				for (int k=0; k<sendLength; ++k) {
					(*this)(sendIndex(k+1)) += u_recv[sendLength-k-1];
				}

				delete[] u_send;
				delete[] u_recv;
				
			}
		}
  	}    

}

template <typename VTYPE>
double*
FLENSDataVector<VTYPE>::vec2c()
{

	double *cVec = new double[(*this).length()];
	
	//Copy FLENS vector to C array:
	for (int i=1; i<=(*this).length(); ++i) {
	
		cVec[i-1] = (*this)(i);
	
	}

	return cVec;

}

}	//namespace flens

namespace flens{ namespace blas{

//Overloaded copy, so that when copying typeII->I, we also apply the appropriate type conversion:
void
copy(FLENSDataVector<FLvTypeII> &orig, FLENSDataVector<FLvTypeI> &dest) 
{
	
	//Create pointers that upcast FLENSDataVector to Parent DenseVector:
	//DenseVector<Array<double> > *tmpOrig = &orig;
	//DenseVector<Array<double> > *tmpDest = &dest;
	
	//Copy data as usual (masquerading as a DenseVector :) ):
	//blas::copy(*tmpOrig, *tmpDest);
	blas::copy(*static_cast<DenseVector<Array<double> > *>(&orig), *static_cast<DenseVector<Array<double> > *>(&dest));

	//Transfer vector type:
	//dest.vType = orig.vType;
	dest.typeII_2_I();

	//(coupling can't be transferred)
}

//Overloaded dot - performs appropriate communication:
double
dot(FLENSDataVector<FLvTypeI> &x1, FLENSDataVector<FLvTypeII> &x2)
{

	//Upcast to DenseVector, and use the standard blas::dot:
	//DenseVector<Array<double> > *tmpx1 = &x1;
	//DenseVector<Array<double> > *tmpx2 = &x2;
		
	double value = blas::dot(*static_cast<DenseVector<Array<double> > *>(&x1), *static_cast<DenseVector<Array<double> > *>(&x2));
	
	//If no communication is required, then we are done:
	if (x1.vType==nonMPI && x2.vType==nonMPI) {
		return value;
	}
	
	//If communication is required..:
	
	//We only multiply typeI with typeII:
	assert(x1.vType != x2.vType);
	
	//Receive buffer:
	double buf = 0;

	//*** Communication to add values from other processes ***/
	MPI::COMM_WORLD.Allreduce(&value, &buf, 1,MPI::DOUBLE,MPI::SUM);
	
	return buf;

}

/*//Overloaded mv, for type updating:
void
mv(Transpose trans, const double &alpha, const GeCRSMatrix<CRS<double, IndexOptions<int, 1> > > &A,
		FLENSDataVector<FlvTypeI> &x, const double &beta, FLENSDataVector<FlvTypeII> &y) {
	
	//Make sure we have the right vector type x:
	assert(x.vType==nonMPI || x.vType==typeI);
	
	//Update y type based on x type:
	if (x.vType == nonMPI) {
		y.vType = nonMPI;
	}
	else {
		y.vType = typeII;
	}
	
	//Use standard blas::mv:
	blas::mv(NoTrans, alpha, A, *static_cast<DenseVector<Array<double> > *>(&x), beta, *static_cast<DenseVector<Array<double> > *>(&y));
		
}*/


}	//namespace blas
}	//namespace flens


#endif	//FLENS_DATA_VECTOR_CPP