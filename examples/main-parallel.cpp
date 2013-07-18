#include <iostream>
#include <iomanip>
#include <mpi.h>

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
    MPI::Init(argc,argv);

    const int rank = MPI::COMM_WORLD.Get_rank();

    if (argc!=2) {
        if (rank==0) {
            std::cerr << "usage: " << argv[0] << "  example_directory" << std::endl;
        }
        MPI::COMM_WORLD.Abort(-1);
    }
    
    IMatrix elements, skeleton;
    IMatrix dirichlet, neumann;
    GeMatrix coordinates;
    IVector elements2procs, sizes(6);
    int numCrossPoints;

    if (rank==0) {
        /* *** only process with rank 0 loads (and sorts) mesh */
        readMeshMPI(argv[1], coordinates, elements, dirichlet, neumann, elements2procs, skeleton, numCrossPoints);
        /* *** first the matrix/vector sizes are sent                */
        sizes(1) =  coordinates.numRows();
        sizes(2) =  elements.numRows();
        sizes(3) =  dirichlet.numRows();
        sizes(4) =  neumann.numRows();
        sizes(5) =  elements2procs.length();
        sizes(6) =  skeleton.numRows();

         MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);
         
    } else {
        MPI::COMM_WORLD.Bcast(sizes.data(), 6, MPI::INT, 0);

           /* *** all other processors allocate their matrices/vectors  */
           coordinates.resize(sizes(1), 2);
           elements.resize(sizes(2), 3);
           dirichlet.resize(sizes(3),2);
           neumann.resize(sizes(4),2);
           elements2procs.resize(sizes(5));
           skeleton.resize(sizes(6), 5);
    }

    /* *** the geometry is distributed  */
    MPI::COMM_WORLD.Bcast(coordinates.data()    , sizes(1)*2, MPI_DOUBLE , 0);
    MPI::COMM_WORLD.Bcast(elements.data()       , sizes(2)*3, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(dirichlet.data()      , sizes(3)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(neumann.data()        , sizes(4)*2, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(elements2procs.data() , sizes(5)  , MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(skeleton.data()       , sizes(6)*5, MPI_INT    , 0);
    MPI::COMM_WORLD.Bcast(&numCrossPoints       , 1         , MPI_INT    , 0);
    
    /* *** create local mesh */
    Mesh mesh(coordinates, elements, dirichlet,neumann, elements2procs, skeleton, numCrossPoints);
    mesh.refineRed();
   	mesh.refineRed();

    mesh.writeData(rank);
    
    
    /* *** create fem object and assemble linear system*/
    FEM<flens::MethMPI> fem(mesh, f, DirichletData,NeumannData); 

    fem.assemble();
    
  fem.solve(cg);
    //fem.solve(gs);

    
	fem.writeSolution(rank);
    
    MPI::Finalize();
    return 0;
}
