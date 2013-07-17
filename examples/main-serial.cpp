#include <iostream>
#include <time.h>
#include <sys/time.h>

#include <flens/flens.cxx>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "../Flens_supl/FlensHeader.h"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"

using namespace std;

typedef flens::GeMatrix<flens::FullStorage<double> > GeMatrix;
typedef flens::GeMatrix<flens::FullStorage<int> > IMatrix;
typedef flens::DenseVector<flens::Array<int> > IVector;


int main(int argc, char *argv[]){
    if (argc!=2) {
        cerr << "usage: " << argv[0] << "  example_directory" << endl;
        exit(-1);
    }
    char directory [256]  = "";
    strcat(directory,argv[1]);

    /*** read the input data */
	IMatrix elements, dirichlet, neumann;
	GeMatrix coordinates;
	cout << "Reading mesh" << endl;
	readMesh(argv[1],coordinates, elements, dirichlet, neumann);
	
	/* *** create mesh */
	cout << "Consturcting mesh" << endl;
	Mesh mesh(coordinates, elements, dirichlet,neumann);
	cout << "Refining" << endl;
	mesh.refineRed();
	mesh.refineRed();
  
  	cout << "Survived mesh" << endl;
    /* *** create fem object and assemble linear system */
    FEM fem(mesh, f, DirichletData,NeumannData);
    cout << "Survived fem" << endl;

    fem.assemble();
    cout << "Survived assemble" << endl;
  
    /* *** write Galerkin matrix */
    //CRSMatrix A = fem.getA(); <------------- Segmentation Fault!!!
    //A.writeFull("./output/A_serial.txt");

	DataVector b = fem.getb();
	b.writeData(0,"./output/b_serial");
	
	/* *** solve problem using the cg method */
    fem.solve(cg);
	
	fem.writeSolution();
	
	
	return 0;
}
