#ifndef MESH_H_
#define MESH_H_

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <mpi.h>
//#include <time.h>

#include <flens/flens.cxx> 
#include "../Flens_supl/CRS_sort.h"
#include "../Flens_supl/FLENS_read_write_methods.h"

#include "Coupling.hpp"


typedef flens::GeMatrix<flens::FullStorage<double> >	GeMatrix;
typedef flens::GeMatrix<flens::FullStorage<int> >     	IMatrix;
typedef flens::DenseVector<flens::Array<double> >     	DVector;
typedef flens::DenseVector<flens::Array<int> >        	IVector;
	
class Mesh{
    
public:

    /* *** public Variables ********************************************************************************/
    GeMatrix coordinates;
    IMatrix elements;
    IMatrix edge2nodes,element2edges;
    std::vector<IVector> dirichlet, neumann, dirichlet2edges, neumann2edges;
    int numEdges, numNodes, numElements;
    int  numDirichlet, numNeumann; // number of connectoed pieces of Dirichlet/Neumann boundary
    
    // variables for MPI
    Coupling coupling;
    
    // variables for adaptive algorithms 
    IVector markedElements;

    
    /* *** constructors and destructor *********************************************************************/
    Mesh();
    Mesh(Mesh &rhs);

    Mesh(const GeMatrix &_coordinates, const IMatrix &_elements,
         IMatrix &_dirichlet, IMatrix &_neumann,
         const IVector _elements2procs = IVector(), 
         const IMatrix _skeleton = IMatrix(), int _numCrossPoints_gl = 0);
       
    ~Mesh();
    
    /* *** operators  *************************************************/
    Mesh & operator=(const Mesh &rhs);
    
    /* *** write mesh for plot *****************************************************************************/                        
    void writeData(int proc=0,std::string dir = "./output/");

    /* *** routines for refinement, (uniform and adaptive) ***************************************************/
    void refineRed();
    IVector refineRGB(DVector &material);

    /* *** check if mesh is distributed (MPI) or not (serial) */
    bool distributed();

    /* *** routines to create all datastructures *************************************************************/
    void provideGeometricData();
    void provideGeometricDataMPI();

private:
    void _distributeMesh(const GeMatrix &coordinates_gl, const IMatrix &elements_gl, 
                       		IMatrix &dirichlet_gl, IMatrix &neumann_gl,
                           	const IVector &elements2procs, const IMatrix &skeleton, 
                           	int numCrossPoints_gl);
                           
    void _convertBdryData(IMatrix bdry_tmp, std::vector<IVector> &bdry);
};
    
/* *** Functions to read the mesh from .dat-files  (serial and parallel) ************************************/
void readMesh(char * input, GeMatrix &coordinates,IMatrix &elements, IMatrix & dirichlet,
              IMatrix &neumann);
              
void readMeshMPI(char * input, GeMatrix &coordinates,IMatrix &elements, IMatrix & dirichlet,
                IMatrix &neumann,IVector &elements2procs, IMatrix &skeleton, int &numCrossPoints);





#endif

