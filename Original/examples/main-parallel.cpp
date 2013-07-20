#include <iostream>
#include <iomanip>
#include <mpi.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"

using namespace std;


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
    
    /* *** Initialise timer object */
    /* *** initialise geometry data containers */
    IndexMatrix elements, skeleton;
    IndexMatrix dirichlet, neumann;
    Matrix coordinates;
    IndexVector elements2procs, sizes(6);
    int numCrossPoints;

	/* *** load input geometry */
    if (rank==0) {
        /* *** only process with rank 0 loads (and sorts) mesh */
        readMeshMPI(argv[1], coordinates, elements, dirichlet, neumann, elements2procs, skeleton, numCrossPoints);
        /* *** first the matrix/vector sizes are sent                */
        sizes(0) =  coordinates.numRows();
        sizes(1) =  elements.numRows();
        sizes(2) =  dirichlet.numRows();
        sizes(3) =  neumann.numRows();
        sizes(4) =  elements2procs.length();
        sizes(5) =  skeleton.numRows();

         MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);
         
    } else {
        MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);

           /* *** all other processors allocate their matrices/vectors  */
           coordinates.resize(sizes(0), 2);
           elements.resize(sizes(1), 3);
           dirichlet.resize(sizes(2),2);
           neumann.resize(sizes(3),2);
           elements2procs.resize(sizes(4));
           skeleton.resize(sizes(5), 5);
    }

    /* *** the geometry is distributed */
    MPI::COMM_WORLD.Bcast(coordinates.data()    , sizes(0)*2, MPI_DOUBLE , 0);
    MPI::COMM_WORLD.Bcast(elements.data()       , sizes(1)*3, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(dirichlet.data()      , sizes(2)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(neumann.data()        , sizes(3)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(elements2procs.data() , sizes(4)  , MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(skeleton.data()       , sizes(5)*5, MPI_INT    , 0);
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
    FEM fem(mesh, f, DirichletData,NeumannData); 

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
