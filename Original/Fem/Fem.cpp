#include "Fem.hpp"

#include "Solver.hpp"

/************************************************************************************************************/
/************************************    constructor   ******************************************************/
/************************************************************************************************************/

/* *** constructor that works for MPI and non-MPI */
FEM::FEM(Mesh &mesh, double (*f)(double,double), 
         double (*DirichletData)(double,double), double (*g)(double,double)):
        _mesh(mesh),
         _A(),
         _uD(mesh.coupling,mesh.numNodes,nonMPI),
         _u(mesh.coupling,mesh.numNodes,nonMPI),
         _b(mesh.coupling,mesh.numNodes,nonMPI)
{
    
    /* *** parallel version             *** */
    /* *** adjust types of data vectors *** */
    if (_mesh.distributed()){    
        _uD.type = typeI;
        _u.type  = typeI;
        _b.type  = typeII;
    }

    /* *** set function pointers */
     _f=f;
     _g=g;
     _DirichletData=DirichletData;
}

FEM::~FEM()
{
    
}


/**********************************************************************************************/
/*******************************    mehtods for assembling   **********************************/
/**********************************************************************************************/
// assemble the Galerkin matrix and the right hand side 
void FEM::assemble()
{
  // write values of dirichlet data 
    _updateDirichlet();
    IndexVector I(9*_mesh.numElements),J(9*_mesh.numElements);
    Vector val(9*_mesh.numElements);
    double c1[2],d21[2],d31[2];
    double area4, fval, a, b, c;
    
    for(int i=0;i<_mesh.numElements;i++){
  
        /* *** assemble stiffness matrix A **********************************************************/
        // compute area of element: 
        // c1 = first vertex of element, d21 vector from first to second vertex, d31 vector 
        // from third to first vertex 
        c1[0]  = _mesh.coordinates( _mesh.elements(i,0)-1,0 );
        c1[1]  = _mesh.coordinates( _mesh.elements(i,0)-1,1 );
        d21[0] = _mesh.coordinates( _mesh.elements(i,1)-1,0 ) - c1[0];
        d21[1] = _mesh.coordinates( _mesh.elements(i,1)-1,1 ) - c1[1];
        d31[0] = _mesh.coordinates( _mesh.elements(i,2)-1,0 ) - c1[0];
        d31[1] = _mesh.coordinates( _mesh.elements(i,2)-1,1 ) - c1[1];
        // compute 4*|T_i| via 2-dimensional cross product
        area4=2*(d21[0]*d31[1]-d21[1]*d31[0]);
        
        // set I and J  (row and column index)
        for(int j=0;j<3;j++){
            I(i*9+j  ) = _mesh.elements(i,j)-1;
            I(i*9+j+3) = _mesh.elements(i,j)-1;
            I(i*9+j+6) = _mesh.elements(i,j)-1;
        
            J(i*9+3*j  ) = _mesh.elements(i,j)-1;
            J(i*9+3*j+1) = _mesh.elements(i,j)-1;
            J(i*9+3*j+2) = _mesh.elements(i,j)-1;
        }
        // set values for stiffness matrix A
        a = (d21[0]*d31[0]+d21[1]*d31[1])/area4;
        b = (d31[0]*d31[0]+d31[1]*d31[1])/area4;
        c = (d21[0]*d21[0]+d21[1]*d21[1])/area4;
    
        val(i*9  ) = -2*a+b+c;
        val(i*9+1) = a-b;
        val(i*9+2) = a-c;
        val(i*9+3) = a-b;
        val(i*9+4) = b;
        val(i*9+5) = -a;
        val(i*9+6) = a-c;
        val(i*9+7) = -a;
        val(i*9+8) = c;
        
        /* *** assemble right-hand side *****************************************************/
        // evaluate volume force f at center of mass of each element
        fval = _f(c1[0]+(d21[0]+d31[0])/3., c1[1]+(d21[1]+d31[1])/3.) /12.*area4;
        _b.values( _mesh.elements(i,0)-1) += fval;
        _b.values( _mesh.elements(i,1)-1) += fval;
        _b.values( _mesh.elements(i,2)-1) += fval;
    }

    
    /* *** incorprate Neumann data */
    /* *** compute bj += int _{Gamma_N} g(x) phi_j(x) dx */
    for(int k=0;k<_mesh.numNeumann;k++){
        //walk through all neumann edges
        for(int j=0;j<_mesh.neumann[k].length()-1;j++){
            double cn1[2],cn2[2], length2;
            cn1[0] = _mesh.coordinates(_mesh.neumann[k](j)  -1,0);
            cn1[1] = _mesh.coordinates(_mesh.neumann[k](j)  -1,1);
            cn2[0] = _mesh.coordinates(_mesh.neumann[k](j+1)-1,0);
            cn2[1] = _mesh.coordinates(_mesh.neumann[k](j+1)-1,1);
            length2=sqrt( (cn1[0]-cn2[0])*(cn1[0]-cn2[0]) + (cn1[1]-cn2[1])*(cn1[1]-cn2[1]) )/2.;
            
            // evaluate Neumann data g at midpoint of neumann edge and multiply with half length 
            double gmE = length2*_g(0.5*(cn1[0]+cn2[0]),0.5*(cn1[1]+cn2[1]));
            _b.values( _mesh.neumann[k](j)  -1) += gmE;
            _b.values( _mesh.neumann[k](j+1)-1) += gmE;    
        }    
    }
    
    /* *** set stiffness matrix */
    _A = CRSMatrix(I,J,val);
    
    /* *** set right-hand side vector b */    
    DataVector Au(_mesh.coupling,_A.numRows(), _b.type);
        
    /* *** incorporate Dirichlet-Data */
    CRSmatVec(Au,_A,_uD);
    add(_b, Au,-1.);    
}

/***********************************************************************************************************/
/**********************************    methods for solving SLE   *******************************************/
/***********************************************************************************************************/
void FEM::solve(Solver method)
{    
    /* *** set fixed nodes (all Dirichletnodes) ***************************************/
    //count number of fixed Nodes
    int numfixed=0;
    for(int k=0;k<_mesh.numDirichlet;k++){
        numfixed += _mesh.dirichlet[k].length();
    } 
    for(int k=0;k<_mesh.coupling.crossPointsBdryData.length();k++){
        if(_mesh.coupling.crossPointsBdryData(k)!=0) numfixed++;        
    }
    
    //set fixedNodes
    IndexVector fixedNodes(numfixed);
    int index=0;
    for(int k=0;k<_mesh.numDirichlet;k++){ 
        fixedNodes.set(index, _mesh.dirichlet[k].length(), _mesh.dirichlet[k].data() );
        index += _mesh.dirichlet[k].length();
    }
    for(int k=0;k<_mesh.coupling.crossPointsBdryData.length();k++){
        if(_mesh.coupling.crossPointsBdryData(k)!=0){
            fixedNodes(index) = k+1;
            index++;
        }
    }
    
    /* *** set maxIt and tolerance *****************************************************/
    maxIt = _mesh.numNodes;
    if(_mesh.coupling.numCoupling>0){
        int maxItgl;
        MPI::COMM_WORLD.Allreduce(&maxIt, &maxItgl, 1, MPI::INT,MPI::SUM);
        maxIt=maxItgl;
    }
    tol = 1e-5;
    
    /* *** solve SLE *********************************************************************/
    int it=0;
    if(method == cg){
        if(!_mesh.distributed()){
            it = CG(_A ,_u.values, _b.values, fixedNodes, maxIt, tol);
        }else{    
        	const int rank = MPI::COMM_WORLD.Get_rank();
        	_b.writeData(rank,"b.txt");
        	if (rank==0){
        	_A.writeCRS("A.txt");}
            it = CG_MPI(_A ,_u, _b, fixedNodes,  maxIt, tol);
        }
    }
    else if (method == gs){
		if(!_mesh.distributed()){
			it = forwardGS( _A, _u.values, _b.values, fixedNodes, maxIt, tol);
		}else{
			it = forwardGS_MPI(_A, _u, _b, _mesh.coupling, fixedNodes, maxIt);
		}
		
	}
	else{
        std::cerr << "Error: Solver method not implemented yet!" << std::endl;
        exit(1);
    }    

    if(!_mesh.distributed() || MPI::COMM_WORLD.Get_rank()==0)
        std::cout <<std::endl<< "Iterations: "<<it<< "  , MaxIterations:  "<<maxIt<<std::endl;

    /* *** update Dirichlet data (x is set to zero at dirichlet nodes in solving methods) ***/
    _updateDirichlet();

}


/**********************************************************************************************************/
/*******************************    methods for refinement   **********************************************/
/**********************************************************************************************************/
void FEM::refineRed()
{
  /* In future projects this function will be modified such that all meshes
     are stored in a hierarchy object which is needed for the multigrid algorithm.
     For now we don't need this and the function just refines the mesh and
     adjusts the size of the vectors 
  */
  
    /* *** refine mesh */
  _mesh.refineRed();

    /* ***resize DataVectors */
  _u.resize(_mesh.numNodes);
  _b.resize(_mesh.numNodes);
  _uD.resize(_mesh.numNodes);
        
}

/**********************************************************************************************/
/*******************************  getter and write methods   **********************************/
/**********************************************************************************************/
CRSMatrix FEM::getA()
{
  return _A;
}

DataVector FEM::getb()
{
    return _b;
}

int FEM::getNumElements()
{
    return _mesh.numElements;
}

void FEM::writeSolution(int proc, std::string filename)
{
  _mesh.writeData(proc, filename);
  _u.writeData(proc, filename + "solution");
}

/**********************************************************************************************/
/*********************************    private methods   ***************************************/
/**********************************************************************************************/

/* write the dirichlet data into the vector _uD and in the solution _u */
void FEM::_updateDirichlet()
{
    int index1;
    for(int k=0; k<_mesh.numDirichlet; k++) {
        for(int i=0;i<_mesh.dirichlet[k].length();i++){
         index1 = _mesh.dirichlet[k](i);
            /* *** set Dirichlet data on nodes */
         _uD.values(index1-1) =  _DirichletData( _mesh.coordinates(index1-1,0), _mesh.coordinates(index1-1,1) );
       /* *** set Dirichlet data in solution vector */
       _u.values(index1-1) = _uD.values(index1-1);
        }
    }
    for(int i=0;i<_mesh.coupling.crossPointsBdryData.length();i++){
        if(_mesh.coupling.crossPointsBdryData(i)!=0){
            /* *** set Dirichlet data on nodes */
      _uD.values(i) = _DirichletData( _mesh.coordinates(i,0), _mesh.coordinates(i,1) );
      /* *** set Dirichlet data in solution vector */
      _u.values(i) = _uD.values(i);
        }
    }    
}        

