#ifndef FEM_H_
#define FEM_H_ 1

#include <flens/flens.cxx>
#include <vector>
#include <math.h>


#include "../Flens_supl/FlensHeader.h"

#include "Mesh.hpp"

#include "../LinearAlgebra/CRSMatrix.hpp"
#include "DataVector.hpp"


// flags for solving SLE
enum Solver { cg, pcg, gs, mg, jac };

template <typename METH>
class FEM{
public:
	
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
	CRSMatrix getA();
	flens::FLENSDataVector<flens::FLvTypeII> getb();
		
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
	flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1> > > fl_A;
	
	flens::FLENSDataVector<typename METH:: I>  fl_uD, fl_u;
	flens::FLENSDataVector<typename METH:: II> fl_b;
		
	//Function pointers for Dirichlet-/Neumann-data and right-hand side:
	double (*_f)(double, double);
	double (*_g)(double, double);
	double (*_DirichletData)(double,double);
		
	/*** Private routines ****************************************************************/
	void _updateDirichlet();
	
};


#include "Fem.cpp"

#endif

