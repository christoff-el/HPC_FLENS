#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <sys/utsname.h>

#include <flens/flens.cxx>

#include "../Flens_supl/FlensHeader.h"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"
#include "timer.h"

using namespace std;

typedef flens::GeMatrix<flens::FullStorage<double> > GeMatrix;
typedef flens::GeMatrix<flens::FullStorage<int> >    IMatrix;
typedef flens::DenseVector<flens::Array<int> >       IVector;


int main(int argc, char *argv[]){

	/* *** initialise MPI interface */
    MPI::Init(argc,argv);
    const int rank = 	   MPI::COMM_WORLD.Get_rank();
    const int procCount =  MPI::COMM_WORLD.Get_size();

	Timer timer;
	timer.start();
	
	/* *** check input parameters */
    if (argc!=5) {
        if (rank==0) {
            std::cerr << "usage: " << argv[0] << "  example_directory, #refinements, solver, timings" << std::endl;
        }
        MPI::COMM_WORLD.Abort(-1);
    }
    
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
    	MPI::COMM_WORLD.Abort(-1);
    	
    }
    
    /* *** output node names, so that we know which nodes/processors/cores of our
    		cluster are actually running the code */
    struct utsname machine_info;
    uname(&machine_info);
    for (int i=0; i<procCount; ++i) {
    	if (rank == i) {
    	
    		MPI::COMM_WORLD.Barrier();
    		cout << "I'm " << machine_info.nodename << 
    				", running the FEM package as rank " << rank << endl;
    				
    	}
    }
    
    /* *** initialise geometry data containers */
    IMatrix elements, skeleton;
    IMatrix dirichlet, neumann;
    GeMatrix coordinates;
    IVector elements2procs, sizes(6);
    int numCrossPoints;

	/* *** load input geometry */
    if (rank==0) {
    
        /* *** only process with rank 0 loads (and sorts) mesh */
        readMeshMPI(argv[1], coordinates, elements, dirichlet, neumann, elements2procs, skeleton, numCrossPoints);
        
        /* *** first the matrix/vector sizes are sent */
        sizes(1) =  coordinates.numRows();
		sizes(2) =  elements.numRows();
        sizes(3) =  dirichlet.numRows();
        sizes(4) =  neumann.numRows();
        sizes(5) =  elements2procs.length();
        sizes(6) =  skeleton.numRows();

		MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);
         
    } 
    else {
		
		/* *** distribute required sizes of data containers from node 0 to world */
		MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);

        /* *** all other processors allocate their matrices/vectors  */
        coordinates.resize(   sizes(1), 2);
        elements.resize(      sizes(2), 3);
        dirichlet.resize(     sizes(3), 2);
    	neumann.resize(       sizes(4), 2);
        elements2procs.resize(sizes(5)   );
		skeleton.resize(      sizes(6), 5);
		
	}

    /* *** the geometry is distributed */
    MPI::COMM_WORLD.Bcast(coordinates.data()    , sizes(1)*2, MPI_DOUBLE , 0);
    MPI::COMM_WORLD.Bcast(elements.data()       , sizes(2)*3, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(dirichlet.data()      , sizes(3)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(neumann.data()        , sizes(4)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(elements2procs.data() , sizes(5)  , MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(skeleton.data()       , sizes(6)*5, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(&numCrossPoints       , 1         , MPI_INT    , 0);
    
    timer.out(rank, "load", timings);
    
    /* *** create local mesh */
    Mesh mesh(coordinates, elements, dirichlet,neumann, elements2procs, skeleton, numCrossPoints);
    
    /* *** refine */
    for (int i=0; i<numRefine; ++i) {
	   	mesh.refineRed();
	}

	timer.out(rank, "mesh", timings);
	
    //mesh.writeData(rank);
    
    /* *** create fem object and assemble linear system*/
    FEM<flens::MethMPI> fem(mesh, f, DirichletData,NeumannData); 

    fem.assemble();
    
    timer.out(rank, "assembly", timings);
    
    /* *** solve problem using specified method */
    fem.solve(solver);
    
    timer.out(rank, "solving", timings);
    
    /* *** each processor writes the solution to its own output file */
	fem.writeSolution(rank);
    
    /* *** finish */
    MPI::Finalize();
    return 0;
    
}
