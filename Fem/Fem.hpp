#ifndef FEM_H_
#define FEM_H_ 1

#include <flens/flens.cxx>
#include <vector>
#include <math.h>
#include <type_traits>


#include "../Flens_supl/FlensHeader.h"

#include "Mesh.hpp"

#include "../LinearAlgebra/CRSMatrix.hpp"
#include "DataVector.hpp"


// flags for solving SLE
enum Solver { cg, pcg, gs, mg, jac };

template <typename A>
struct SelectDataVector
{
};

template <>
struct SelectDataVector<flens::MethMPI>
{
    typedef flens::FLENSDataVector<flens::FLvTypeI> 	TypeI;
    typedef flens::FLENSDataVector<flens::FLvTypeII> 	TypeII;
};

template <>
struct SelectDataVector<flens::MethNonMPI>
{
    typedef flens::FLENSDataVector<flens::FLvNonMPI> 	TypeI;
    typedef flens::FLENSDataVector<flens::FLvNonMPI> 	TypeII;
};


template <typename METH>
class FEM{
public:

	typedef flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1> > >	CRSMat;

	/*** Constructors and destructor *****************************************************/
	FEM(Mesh &mesh_tmp, 
			double (*f)(double,double), 
			double (*DirichletData)(double,double), 
			double (*g)(double,double));
	
	~FEM();	
		
	/*** Method for assembly *************************************************************/
	void assemble();
		
	/*** Method for solving SLE **********************************************************/
	void solve(Solver method);
		
	/*** Method for refinement ***********************************************************/
	void refineRed();
		
	/*** getter and write methods ********************************************************/
	CRSMat getA();
	typename SelectDataVector<METH>::TypeII getb();

		
    int getNumElements();
	void writeSolution(int proc=0,std::string filename="./output/");		
			
	/*** Public variables ****************************************************************/
	
	//Parameters for solving SLE:
	int maxIt;
	double tol; 
		
private:

	/*** Private variables ***************************************************************/
	
	//Local mesh:
	Mesh _mesh; 
		
	//Storage structures for FEM system:
	CRSMat 	_A;

	
	typename SelectDataVector<METH>::TypeI    _uD, _u;
	typename SelectDataVector<METH>::TypeII   _b;
		
	//Function pointers for Dirichlet-/Neumann-data and right-hand side:
	double (*_f)(double, double);
	double (*_g)(double, double);
	double (*_DirichletData)(double,double);
		
	/*** Private routines ****************************************************************/
	void _updateDirichlet();
	
};


#include "Fem.cpp"

#endif

