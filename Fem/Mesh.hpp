#ifndef MESH_H_
#define MESH_H_

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <mpi.h>
#include <time.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "Coupling.hpp"

class Mesh{
    
public:
    /* *** public Variables ********************************************************************************/
    Matrix coordinates;
    IndexMatrix elements;
    IndexMatrix edge2nodes,element2edges;
    std::vector<IndexVector> dirichlet, neumann, dirichlet2edges, neumann2edges;
    int numEdges, numNodes, numElements;
    int  numDirichlet, numNeumann; // number of connectoed pieces of Dirichlet/Neumann boundary
    // variables for MPI
    Coupling coupling;
    // variables for adaptive algorithms 
    IndexVector markedElements;

    
    /* *** constructors and destructor *********************************************************************/
    Mesh();
    Mesh(Mesh &rhs);

    Mesh(const Matrix &_coordinates, const IndexMatrix &_elements,
         const IndexMatrix &_dirichlet, const IndexMatrix &_neumann,
         const IndexVector _elements2procs = IndexVector(), const IndexMatrix _skeleton = IndexMatrix(), int _numCrossPoints_gl = 0);
       
    ~Mesh();
    
    /* *** operators  *************************************************/
    Mesh & operator=(const Mesh &rhs);
    
    /* *** write mesh for plot *******************************************************************************/                        
    void writeData(int proc=0,std::string dir = "./output/");

    /* *** routines for refinement, (uniform and adaptive) ***************************************************/
    void refineRed();
    IndexVector refineRGB(Vector &material);

    /* *** check if mesh is distributed (MPI) or not (serial) */
    bool distributed();

    /* *** routines to create all datastructures *************************************************************/
    void provideGeometricData();
    void provideGeometricDataMPI();

private:
    void _distributeMesh(const Matrix &coordinates_gl, const IndexMatrix &elements_gl, 
                       const IndexMatrix &dirichlet_gl, const IndexMatrix &neumann_gl,
                           const IndexVector &elements2procs, const IndexMatrix &skeleton, int numCrossPoints_gl);
                           
    void _convertBdryData(const IndexMatrix bdry_tmp, std::vector<IndexVector> &bdry);
};
    
/* *** Functions to read the mesh from .dat-files  (serial and parallel) ************************************/
void readMesh(char * input, Matrix &coordinates,IndexMatrix &elements, IndexMatrix & dirichlet,
              IndexMatrix &neumann);
void readMeshMPI(char * input, Matrix &coordinates,IndexMatrix &elements, IndexMatrix & dirichlet,
                IndexMatrix &neumann,IndexVector &elements2procs, IndexMatrix &skeleton, int &numCrossPoints);





#endif

