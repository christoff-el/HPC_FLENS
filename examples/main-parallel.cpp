#include <iostream>
#include <iomanip>
#include <mpi.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "../Flens_supl/FlensHeader.h"
#include "../Fem/FemHeader.hpp"
#include "functions.hpp"

using namespace std;


int main(int argc, char *argv[]){
    MPI::Init(argc,argv);

    const int rank = MPI::COMM_WORLD.Get_rank();

    if (argc!=2) {
        if (rank==0) {
            std::cerr << "usage: " << argv[0] << "  example_directory" << std::endl;
        }
        MPI::COMM_WORLD.Abort(-1);
    }
    
    IndexMatrix elements, skeleton;
    IndexMatrix dirichlet, neumann;
    Matrix coordinates;
    IndexVector elements2procs, sizes(6);
    int numCrossPoints;
    
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

    /* *** the geometry is distributed  */
    MPI::COMM_WORLD.Bcast(coordinates.data()    , sizes(0)*2, MPI_DOUBLE , 0);
    MPI::COMM_WORLD.Bcast(elements.data()       , sizes(1)*3, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(dirichlet.data()      , sizes(2)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(neumann.data()        , sizes(3)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(elements2procs.data() , sizes(4)  , MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(skeleton.data()       , sizes(5)*5, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(&numCrossPoints       , 1         , MPI_INT    , 0);
    
    /* *** create local mesh */
    Mesh mesh(coordinates, elements, dirichlet,neumann, elements2procs, skeleton, numCrossPoints);
    mesh.refineRed();
    mesh.refineRed();
    mesh.refineRed();
    
    mesh.writeData(rank);
    
    
    /* *** create fem object and assemble linear system*/
    FEM fem(mesh, f, DirichletData,NeumannData); 
    fem.assemble();

    fem.solve(cg);
	fem.writeSolution(rank);
    
    MPI::Finalize();
    return 0;
}
