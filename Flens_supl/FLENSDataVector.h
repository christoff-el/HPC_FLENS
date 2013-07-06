#ifndef FLENS_DATA_VECTOR_H
#define FLENS_DATA_VECTOR_H 1

#include <mpi.h>
#include <flens/flens.cxx>

//#include "flens2c.h"

#include "../Fem/Coupling.hpp"


namespace flens{

enum VectorType {typeI, typeII, nonMPI};

struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	//MPI Version (with Coupling, requires VectorType):
	explicit
    FLENSDataVector(int n, const Coupling &_coupling, const VectorType _vType);
    
    //Non-MPI version (no coupling, assumes 'nonMPI' VectorType):
    explicit
    FLENSDataVector(int n);
    
    //Copy constructor:
    explicit
    FLENSDataVector(const FLENSDataVector &rhs);
    

	//Member objects:
    VectorType vType;
	const Coupling &coupling;

	//Member methods:
	void typeII_2_I();
	void typeI_2_II();
	
	void commCrossPoints();
	void commBoundaryNodes();
	
	double* vec2c();
	

};


}	//namespace flens


#endif	//FLENS_DATA_VECTOR_H