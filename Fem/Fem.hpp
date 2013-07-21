#ifndef FEM_H_
#define FEM_H_ 1

#include <flens/flens.cxx>
#include <math.h>
#include <type_traits>

#include "../Flens_supl/FlensHeader.h"
#include "Mesh.hpp"


// flags for solving SLE
enum Solver { cg, pcg, gs, mg, jac };

//Structure for selecting either MPI or NonMPI DataVectors:
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
	typedef flens::DenseVector<flens::Array<double> >								DenseVector;

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
	CRSMat 		_A;
	DenseVector _uD;
	
	typename SelectDataVector<METH>::TypeI    _u;			//Vector types based on FEM METHOD
	typename SelectDataVector<METH>::TypeII   _b;			//    type (MPI or NonMPI).
		
	//Function pointers for Dirichlet-/Neumann-data and right-hand side:
	double (*_f)(double, double);
	double (*_g)(double, double);
	double (*_DirichletData)(double,double);
		
	/*** Private routines ****************************************************************/
	void _updateDirichlet();
	
};


#include "Fem.cpp"

#endif	//FEM_H

