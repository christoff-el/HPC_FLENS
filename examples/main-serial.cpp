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

	Timer timer;
	timer.start();
	
	/* *** check input parameters */
    if (argc!=5) {
        cerr << "usage: " << argv[0] << "  example_directory, #refinements, solver, timings" << endl;
        exit(-1);
    }
    char directory [256]  = "";
    strcat(directory,argv[1]);
    
    int numRefine = atoi(argv[2]);
    bool timings =  atoi(argv[4]);
    
    /* *** parse solver input */
    Solver solver;
    string chosen_solver = argv[3];
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
	
	timer.out(0, "load", timings);
	
	/* *** create and refine mesh */
	Mesh mesh(coordinates, elements, dirichlet,neumann);
	
	/* *** refine */
    for (int i=0; i<numRefine; ++i) {
	   	mesh.refineRed();
	}
	
	timer.out(rank, "mesh", timings);
  
    /* *** create fem object and assemble linear system */
    FEM<flens::MethNonMPI> fem(mesh, f, DirichletData, NeumannData);

    fem.assemble();
    
    timer.out(rank, "assembly", timings);
	
	/* *** solve problem using specified method */
    fem.solve(solver);
    
    timer.out(rank, "solving", timings);

	/* *** write solution to file */
	fem.writeSolution();
	
	
	return 0;
}
