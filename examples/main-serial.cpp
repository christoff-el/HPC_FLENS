#include <iostream>
#include <time.h>
#include <sys/time.h>

#include <flens/flens.cxx>

#include "../Flens_supl/FlensHeader.h"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"

using namespace std;

typedef flens::GeMatrix<flens::FullStorage<double> > GeMatrix;
typedef flens::GeMatrix<flens::FullStorage<int> > 	 IMatrix;



int main(int argc, char *argv[]){

	/* *** check input parameters */
    if (argc!=3) {
        cerr << "usage: " << argv[0] << "  example_directory, solver" << endl;
        exit(-1);
    }
    char directory [256]  = "";
    strcat(directory,argv[1]);
    
    /* *** parse solver input */
    Solver solver;
    string chosen_solver = argv[2];
     if (chosen_solver == "cg" || chosen_solver == "CG") {
    
    	solver = cg;
    	
    }
    else if (chosen_solver == "gs" || chosen_solver == "GS") {
    
    	solver = gs;
    	
    }
    else {
    
    	cerr << "Invalid solver specified!" << endl;
    	exit(-1);
    	
    }    

    /*** read the input data */
	IMatrix  elements, dirichlet, neumann;
	GeMatrix coordinates;
	readMesh(argv[1],coordinates, elements, dirichlet, neumann);
	
	/* *** create and refine mesh */
	Mesh mesh(coordinates, elements, dirichlet,neumann);
	mesh.refineRed();
	mesh.refineRed();
  
    /* *** create fem object and assemble linear system */
    FEM<flens::MethNonMPI> fem(mesh, f, DirichletData, NeumannData);

    fem.assemble();
	
	/* *** solve problem using specified method */
    fem.solve(solver);

	/* *** write solution to file */
	fem.writeSolution();
	
	
	return 0;
}
