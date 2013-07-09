#ifndef FEM_H_
#define FEM_H_

#include <flens/flens.cxx>
#include <vector>
#include <math.h>

//#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "../LinearAlgebra/CRSMatrix.hpp"
#include "../Flens_supl/FlensHeader.h"
//#include "MathOperationsMPI.hpp"
#include "Mesh.hpp"
#include "DataVector.hpp"

//using namespace flens;


// flags for solving SLE
enum Solver { cg, pcg, gs, mg, jac };

class FEM{
public:
	
	/* *** constructors and destructor *****************************************************/
	FEM(Mesh &mesh_tmp, 
	double (*f)(double,double), 
	double (*DirichletData)(double,double), 
	double (*g)(double,double));
	
	~FEM();	
		
	/* *** methods for assembling **********************************************************/
	void assemble();
		
	/* *** methods for solving SLE *********************************************************/
	void solve(Solver method);
		
	/* *** methods for refinement **********************************************************/
	void refineRed();
		
	/* *** getter and write methods ********************************************************/
	CRSMatrix getA();
	DataVector getb();
		
    int getNumElements();
	void writeSolution(int proc=0,std::string filename="./output/");		
			
	/* *** public variables ****************************************************************/
	// variables for solving SLE
	int maxIt;
	double tol; 
		
private:
	/* *** private variables **************************************************************/
	// local mesh
	Mesh _mesh; 
		
	// variables for FEM
	flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1> > > fl_A;
	flens::FLENSDataVector fl_uD, fl_u, fl_b;
		
	// function pointers for Dirichlet-/Neumann-data and right-hand side
	double (*_f)(double, double);
	double (*_g)(double, double);
	double (*_DirichletData)(double,double);
		
	/* *** private routines **************************************************************/
	void _updateDirichlet();
	
};

#endif

