#ifndef FLENS_DATA_VECTOR_H
#define FLENS_DATA_VECTOR_H 1

#include <mpi.h>
#include <flens/flens.cxx>

#include "../Fem/Coupling.hpp"


namespace flens{

class FLvNonMPI;
class FLvTypeI;
class FLvTypeII;

struct MethMPI {
	typedef FLvTypeI  I;
	typedef FLvTypeII II;
};

struct MethNonMPI {
	typedef FLvNonMPI I;
	typedef FLvNonMPI II;
};

template <typename VTYPE = FLvNonMPI>
struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	//Non-MPI constructor:
	FLENSDataVector(int n);
	
	//MPI constructor:
	FLENSDataVector(int n, const Coupling &_coupling);
    
	//Copy constructor (VTYPEs obligated to match):
	FLENSDataVector(const FLENSDataVector<VTYPE> &rhs);
	    

	//Member objects:
    //VectorType vType;
	const Coupling &coupling;

	//Member methods:
	void typeII_2_I();				// <-- applies conversion within a TypeI	(data already copied).
	void typeI_2_II();				// <-- applies conversion within a TypeII
	
	void commCrossPoints();
	void commBoundaryNodes();
	
	void writeData(int proc, std::string filename);
	double* vec2c();

};

}	//namespace flens

namespace flens{ namespace blas{

void 
copy(FLENSDataVector<FLvTypeII> &orig, FLENSDataVector<FLvTypeI> &dest);

double
dot(FLENSDataVector<FLvTypeI> &x1, FLENSDataVector<FLvTypeII> &x2);

double
dot(FLENSDataVector<FLvTypeII> &x1, FLENSDataVector<FLvTypeI> &x2);


}	//namespace blas
}	//namespace flens


#include "FLENSDataVector.cpp"

#endif	//FLENS_DATA_VECTOR_H