#include <iostream>
#include <time.h>
#include <sys/time.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"

using namespace std;


int main(int argc, char *argv[]){
    if (argc!=2) {
        cerr << "usage: " << argv[0] << "  example_directory" << endl;
        exit(-1);
    }
    char directory [256]  = "";
    strcat(directory,argv[1]);

    /*** read the input data */
	IndexMatrix elements, dirichlet, neumann;
	Matrix coordinates;
	readMesh(argv[1],coordinates, elements, dirichlet, neumann);
	
	/* *** create mesh */
	Mesh mesh(coordinates, elements, dirichlet,neumann);
  
    /* *** create fem object and assemble linear system */
    FEM fem(mesh, f, DirichletData,NeumannData);
    fem.assemble();
  
    /* *** write Galerkin matrix */
    CRSMatrix A = fem.getA();
    A.writeFull("./output/A_serial.txt");
	
	DataVector b = fem.getb();
	b.writeData(0,"./output/b_serial");
	
	/* *** solve problem using the cg method */
    fem.solve(cg);
	
	fem.writeSolution();
	
	
	return 0;
}
