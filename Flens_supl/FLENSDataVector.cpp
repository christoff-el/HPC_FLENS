#ifndef FLENS_DATA_VECTOR_CPP
#define FLENS_DATA_VECTOR_CPP 1

#include "FLENSDataVector.h"

namespace flens{


//Non-MPI --> Restricted to only FLnonMPI specialisation.
template <typename VTYPE>
FLENSDataVector<VTYPE>::FLENSDataVector(int n)
	:	coupling(Coupling())
{
	VTYPE::CHK;			//<-- If scope ever reaches here, compilation will fail.
						//		e.g. if FLENSDataVector<double> instantiated.
}

template <>
FLENSDataVector<FLvNonMPI>::FLENSDataVector(int n)
	: 	DenseVector<Array<double> >(n),
		coupling(Coupling())
{
	//Permits instatiation of FlNonMPI specialisation.
}


//MPI --> Restricted to only FLvTypeI, FLvTypeII specialisations. 
//				FLvNonMPI is permitted, but supplied 
template <typename VTYPE>
FLENSDataVector<VTYPE>::FLENSDataVector(int n, const Coupling &_coupling)
	:	coupling(Coupling())
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


//Copy constructor (VTYPEs obligated to match):
template <typename VTYPE>
FLENSDataVector<VTYPE>::FLENSDataVector(const FLENSDataVector<VTYPE> &rhs)
   	:	DenseVector<Array<double> >(rhs),			//copy data via flens framework
		coupling(rhs.coupling)
{
}

template <>	
void
FLENSDataVector<FLvTypeI>::typeII_2_I()
{

	//Sum up values at cross points:
	commCrossPoints();

	//Sum up values at boundary nodes:
	commBoundaryNodes();

}

template <>
void
FLENSDataVector<FLvTypeI>::typeI_2_II()
{

	//Divide values at cross points by the number of processes:
	for (int i=1; i<=coupling.local2globalCrossPoints.length(); ++i) {
		(*this)(i) /= coupling.crossPointsNumProcs(i);
	}

	//Divide values at boundary nodes by 2 (since here 2 processes overlap):		(don't divide for cross points!)
	for (int i=1; i<=coupling.numCoupling; ++i) {
		for (int j=2; j<=coupling.boundaryNodes[i].length()-1; ++j) {
			(*this)(coupling.boundaryNodes[i](j)) /= 2.;
		}
	}

}

template <typename VTYPE>
void
FLENSDataVector<VTYPE>::commCrossPoints()
{

	DenseVector<Array<double> > u_crossPoints(coupling.numCrossPoints);
	DenseVector<Array<double> > u_crossPoints_gl(coupling.numCrossPoints);

	//Local values at all global cross points:
	for (int i=1; i<=coupling.local2globalCrossPoints.length(); ++i) {
		u_crossPoints(coupling.local2globalCrossPoints(i)) = (*this)(i);
	}

	/*** MPI Communication for global cross points ***/
	MPI::COMM_WORLD.Allreduce(u_crossPoints.data(),
						      u_crossPoints_gl.data(),
							  coupling.numCrossPoints,
							  MPI::DOUBLE, MPI::SUM);

	for (int i=1; i<=coupling.local2globalCrossPoints.length(); ++i) {
		(*this)(i) = u_crossPoints_gl(coupling.local2globalCrossPoints(i));
	}

}

template <typename VTYPE>
void 
FLENSDataVector<VTYPE>::commBoundaryNodes()
{

	for (int i=1; i<=coupling.maxColor; ++i) {
		for (int j=0; j<coupling.numCoupling; ++j) {
			if (coupling.colors(j+1) == i && coupling.boundaryNodes[j].length()-2 > 0) {

				//Only communicate if there is a boundary node on coupling boundary (no cross points):
				int sendLength = coupling.boundaryNodes[j].length()-2;

				DenseVector<Array<double> > u_send(sendLength);
				DenseVector<Array<double> > u_recv(sendLength);

				//Set local values to be sent:
				for (int k=1; k<=sendLength; ++k) {
					u_send(k) = (*this)(coupling.boundaryNodes[j](k+1));
				}

				//Get values from other processes:
				MPI::COMM_WORLD.Sendrecv(u_send.data(), sendLength, MPI::DOUBLE, 
											coupling.neighbourProcs(j+1)-1, 0,
											u_recv.data(), sendLength, MPI::DOUBLE, 
											coupling.neighbourProcs(j+1)-1,0);

				//Add values collected from the other processes (!!numbering is opposite):
				for (int k=0; k<sendLength; ++k) {
					(*this)(coupling.boundaryNodes[j](k+2)) += u_recv(sendLength-k);
				}
			}
		}
  	}    

}

template <typename VTYPE>
void 
FLENSDataVector<VTYPE>::writeData(int proc, std::string filename)
{

    std::string strproc;
  
    if (proc==0)  strproc="";
    else {
    	std::stringstream ss;
    	ss << proc;
    	strproc = ss.str();
    }
    
    filename = filename + strproc + ".dat";

	std::fstream f;

	f.open(filename.c_str(), std::ios::out);
	if (f.is_open()){

		for (int i=1; i<=(*this).length(); ++i) {
			f << (*this)(i) << std::endl;
		}

		f.close();

	}
    
}

}	//namespace flens



namespace flens{ namespace blas{

//Overloaded copy, so that when copying typeII->I, we also apply the appropriate type conversion:
void
copy(FLENSDataVector<FLvTypeII> &orig, FLENSDataVector<FLvTypeI> &dest) 
{

	//Copy data as usual (masquerading as a DenseVector :) ):
	blas::copy(*static_cast<DenseVector<Array<double> > *>(&orig),
			   *static_cast<DenseVector<Array<double> > *>(&dest));

	//Perform vector type conversion:
	dest.typeII_2_I();
	//(coupling can't be transferred)
}

//Overloaded dot - performs appropriate communication:
double
dot(FLENSDataVector<FLvTypeI> &x1, FLENSDataVector<FLvTypeII> &x2)
{

	//Upcast to DenseVector, and use the standard blas::dot:
	double value = blas::dot(*static_cast<DenseVector<Array<double> > *>(&x1),
	                         *static_cast<DenseVector<Array<double> > *>(&x2));

	//Receive buffer:
	double buf = 0;

	//*** Communication to add values from other processes ***/
	MPI::COMM_WORLD.Allreduce(&value, &buf, 1,MPI::DOUBLE,MPI::SUM);

	return buf;

}

//Adds commutativity to dot:
double
dot(FLENSDataVector<FLvTypeII> &x1, FLENSDataVector<FLvTypeI> &x2)
{
	
	return dot(x2,x1);
	
}


}	//namespace blas
}	//namespace flens


#endif	//FLENS_DATA_VECTOR_CPP