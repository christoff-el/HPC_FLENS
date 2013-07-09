#include "Fem.hpp"

#include "Solver.hpp"


/************************************************************************************************************/
/************************************    constructor   ******************************************************/
/************************************************************************************************************/

/* *** constructor that works for MPI and non-MPI */
FEM::FEM(Mesh &mesh, double (*f)(double,double), 
        		double (*DirichletData)(double,double), double (*g)(double,double))
        :	_mesh(mesh),
        	fl_A(),
	    	fl_uD(mesh.numNodes, mesh.coupling, (flens::VectorType)nonMPI),
         	fl_u(mesh.numNodes, mesh.coupling, (flens::VectorType)nonMPI),
         	fl_b(mesh.numNodes, mesh.coupling, (flens::VectorType)nonMPI),
         	_f(f),
         	_g(g),
         	_DirichletData(DirichletData)
{
    
    //^^^ We initialise an empty CRS matrix via an empty coordinate storage sparse matrix.
    
    /*** Parallel version ***/
    if (_mesh.distributed()) {  

        fl_uD.vType = (flens::VectorType)typeI;
        fl_u.vType  = (flens::VectorType)typeI;
        fl_b.vType  = (flens::VectorType)typeII;
    }

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
  	//Write values of dirichlet data:
    _updateDirichlet();
    
    flens::DenseVector<flens::Array<int> > 	I(9*_mesh.numElements);
    flens::DenseVector<flens::Array<int> > 	J(9*_mesh.numElements);
    flens::DenseVector<flens::Array<double> > val(9*_mesh.numElements);
    double c1[2], d21[2], d31[2];
    double area4, fval, a, b, c;
    
    for (int i=0; i<_mesh.numElements; ++i) {
  
        /*** Assemble stiffness matrix A ***/
        
        //Compute area of element: 
        // c1 = first vertex of element, 
        // d21 vector from first to second vertex, 
        // d31 vector from third to first vertex:
        c1[0]  = _mesh.coordinates( _mesh.elements(i,0)-1,0 );
        c1[1]  = _mesh.coordinates( _mesh.elements(i,0)-1,1 );
        d21[0] = _mesh.coordinates( _mesh.elements(i,1)-1,0 ) - c1[0];
        d21[1] = _mesh.coordinates( _mesh.elements(i,1)-1,1 ) - c1[1];
        d31[0] = _mesh.coordinates( _mesh.elements(i,2)-1,0 ) - c1[0];
        d31[1] = _mesh.coordinates( _mesh.elements(i,2)-1,1 ) - c1[1];
        
        //Compute 4*|T_i| via 2-dimensional cross product:
        area4 = 2*(d21[0]*d31[1] - d21[1]*d31[0]);
        
        //Set I and J  (row and column index)
        for(int j=0; j<3; ++j) {
        
            I(i*9+j  +1) = _mesh.elements(i,j);
            I(i*9+j+3+1) = _mesh.elements(i,j);
            I(i*9+j+6+1) = _mesh.elements(i,j);
        
            J(i*9+3*j  +1) = _mesh.elements(i,j);
            J(i*9+3*j+1+1) = _mesh.elements(i,j);
            J(i*9+3*j+2+1) = _mesh.elements(i,j);
            
        }
        
        //Set values for stiffness matrix A:
        a = (d21[0]*d31[0] + d21[1]*d31[1]) / area4;
        b = (d31[0]*d31[0] + d31[1]*d31[1]) / area4;
        c = (d21[0]*d21[0] + d21[1]*d21[1]) / area4;
    
        val(i*9  +1) = -2*a + b+c;
        val(i*9+1+1) = a-b;
        val(i*9+2+1) = a-c;
        val(i*9+3+1) = a-b;
        val(i*9+4+1) = b;
        val(i*9+5+1) = -a;
        val(i*9+6+1) = a-c;
        val(i*9+7+1) = -a;
        val(i*9+8+1) = c;
        
        
        /*** Assemble right-hand side ***/
        
        //Evaluate volume force f at center of mass of each element:
        fval = _f( c1[0] + (d21[0]+d31[0])/3., c1[1] + (d21[1]+d31[1])/3. ) / 12. * area4;
        fl_b( _mesh.elements(i,0)) += fval;
        fl_b( _mesh.elements(i,1)) += fval;
        fl_b( _mesh.elements(i,2)) += fval;
        
    }

    
    /*** Incorprate Neumann data ***/

    //Compute bj += int _{Gamma_N} g(x) phi_j(x) dx:
    for(int k=0; k<_mesh.numNeumann; ++k) {
    
        //Walk through all neumann edges:
        for(int j=0; j<_mesh.neumann[k].length()-1; ++j) {
        
            double cn1[2],cn2[2], length2;
            cn1[0]  = _mesh.coordinates(_mesh.neumann[k](j)  -1,0);
            cn1[1]  = _mesh.coordinates(_mesh.neumann[k](j)  -1,1);
            cn2[0]  = _mesh.coordinates(_mesh.neumann[k](j+1)-1,0);
            cn2[1]  = _mesh.coordinates(_mesh.neumann[k](j+1)-1,1);
            length2 = sqrt( (cn1[0]-cn2[0])*(cn1[0]-cn2[0]) + (cn1[1]-cn2[1])*(cn1[1]-cn2[1]) ) / 2.;
            
            //Evaluate Neumann data g at midpoint of neumann edge and multiply with half length:
            double gmE = length2*_g(0.5*(cn1[0]+cn2[0]) , 0.5*(cn1[1]+cn2[1]));
            fl_b( _mesh.neumann[k](j)  ) += gmE;
            fl_b( _mesh.neumann[k](j+1)) += gmE;   
             
        }   
         
    }


    /*** Set stiffness matrix ***/
    
    //Build FLENS CRS matrix from I, J, vals (rows, cols, values):
    flens::GeCoordMatrix<flens::CoordStorage<double> > fl_A_coord(fl_uD.length(),fl_uD.length());
    //, flens::CoordRowColCmp, flens::IndexBaseOne<int> 
    for (int i=1; i<=I.length(); ++i) {
    
    	fl_A_coord(I(i),J(i)) += val(i);
    	
    }
    
    fl_A = fl_A_coord;


    /*** Set right-hand side vector b ***/    
    flens::FLENSDataVector fl_Au(fl_uD.length(), _mesh.coupling, fl_b.vType);
        

    /*** Incorporate Dirichlet-Data ***/
    flens::blas::mv(flens::NoTrans, 1., fl_A, fl_uD, 0., fl_Au);
    flens::blas::axpy(-1., fl_Au, fl_b);    

}


/***********************************************************************************************************/
/**********************************    methods for solving SLE   *******************************************/
/***********************************************************************************************************/
void FEM::solve(Solver method)
{    

    /*** Set fixed nodes (all Dirichletnodes) ***/
    
    //Count the number of fixed nodes:
    int numfixed = 0;
    for (int k=0; k<_mesh.numDirichlet; ++k) {
        numfixed += _mesh.dirichlet[k].length();
    } 
    
    for (int k=0; k<_mesh.coupling.crossPointsBdryData.length(); ++k) {
        if(_mesh.coupling.crossPointsBdryData(k)!=0) numfixed++;        
    }
    
    //Set fixedNodes:
    flens::DenseVector<flens::Array<int, flens::IndexOptions<int, 0> > > fixedNodes(numfixed);		//!!Base 0
    int index=0;
    for (int k=0; k<_mesh.numDirichlet; ++k) { 
    
    	//Setting:
    	for (int l=index; l<index+_mesh.dirichlet[k].length(); ++l) {
    		fixedNodes(l) = _mesh.dirichlet[k](l);			//fixedNodes.set(index, _mesh.dirichlet[k].length(), _mesh.dirichlet[k].data() );
        }
        
        index += _mesh.dirichlet[k].length();
        
    }
    
    for (int k=0; k<_mesh.coupling.crossPointsBdryData.length(); ++k) {
    
        if (_mesh.coupling.crossPointsBdryData(k) != 0) {
        
            fixedNodes(index) = k+1;
            index++;
            
        }
        
    }
    
    
    /*** Set maxIt and tolerance ***/
    
    //Set maximum iteration count to the number of nodes:
    maxIt = _mesh.numNodes;
    
    //Distribute maxIt if using MPI:
    if (_mesh.coupling.numCoupling>0) {
    
        int maxItgl;
        MPI::COMM_WORLD.Allreduce(&maxIt, &maxItgl, 1, MPI::INT,MPI::SUM);
        maxIt=maxItgl;
        
    }
    
    //Set tolerance:
    tol = 1e-5;
    
    
    /*** Solve the SLE ***/
    
    int it=0;
    
    //CG Solver:
    if (method == cg) {
    	
    	//Serial solver:
        if (!_mesh.distributed()) {
            it = cg_nompi_blas(fl_A ,fl_u, fl_b, fixedNodes, maxIt, tol);
        }
        
        //Parallel solver:
        else { 
        	const int rank = MPI::COMM_WORLD.Get_rank();  
        	if(rank==0){
        	std::cout<<fl_b<<std::endl;
        	std::cout<<std::endl;
        	std::cout<<fl_A<<std::endl;}
            it = cg_mpi_blas(fl_A ,fl_u, fl_b, fixedNodes,  maxIt, tol);
        }
        
    }
    
    //GS Solver:
    else if (method == gs) {
    
    	//HACK until GS implemented for flens:
    	//Convert u, b back to Funken:			TURNED OFF, since requires a FLENS -> CRS converter
    	
    	//DataVector _u(fl_u.coupling, fl_u.length(), fl_u.vType)
    	//DataVector _b(fl_b.coupling, fl_b.length(), fl_b.vType);
    	//flens2funk_DataVector(fl_u, _u);
    	//flens2funk_DataVector(fl_b, _b);
    	
    	//Serial solver:
		if (!_mesh.distributed()) {
			it = forwardGS(fl_A, fl_u, fl_b, fixedNodes, maxIt, tol);
		}
		
		//Parallel solver:
		else {
			it = forwardGS_MPI(fl_A, fl_u, fl_b, _mesh.coupling, fixedNodes, maxIt);
		}
		
	}
	
	//Something else!:
	else {
	
        std::cerr << "Error: Solver method not implemented yet!" << std::endl;
        exit(1);
        
    }    

	//Only rank 0 (or serial implementation) output required iterations:
    if (!_mesh.distributed() || MPI::COMM_WORLD.Get_rank()==0) {
    
        std::cout <<std::endl<< "Iterations: "<<it<< "  , MaxIterations:  "<<maxIt<<std::endl;

	}
	
	
    /*** Update Dirichlet data (x is set to zero at dirichlet nodes in solving methods) ***/
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
  
    /*** Refine mesh ***/
	_mesh.refineRed();

    /*** Resize DataVectors (note: clears data) ***/
  	fl_u.resize(_mesh.numNodes);
  	fl_b.resize(_mesh.numNodes);
  	fl_uD.resize(_mesh.numNodes);
        
}

/**********************************************************************************************/
/*******************************  getter and write methods   **********************************/
/**********************************************************************************************/
CRSMatrix FEM::getA()
{
	//HACK Err.. Turn if off for now.
  	return CRSMatrix();
}

DataVector FEM::getb()
{
	//HACK to return b as a DataVector:
	DataVector _b(fl_b.coupling, fl_b.length(), (vectorType)fl_b.vType);
	flens2funk_DataVector(fl_b, _b);
	
    return _b;
}

int FEM::getNumElements()
{
    return _mesh.numElements;
}

void FEM::writeSolution(int proc, std::string filename)
{
	_mesh.writeData(proc, filename);
	
	//HACK to write from DataVector:
	DataVector _u(fl_u.coupling, fl_u.length(), (vectorType)fl_u.vType);
	flens2funk_DataVector(fl_u, _u);
	
	_u.writeData(proc, filename + "solution");
	
}

/**********************************************************************************************/
/*********************************    private methods   ***************************************/
/**********************************************************************************************/

/* write the dirichlet data into the vector _uD and in the solution _u */
void FEM::_updateDirichlet()
{

    int index1;
    for (int k=0; k<_mesh.numDirichlet; ++k) {
        for (int i=0; i<_mesh.dirichlet[k].length(); ++i) {
        
        	index1 = _mesh.dirichlet[k](i);
        	
            //Set Dirichlet data on nodes:
         	fl_uD(index1) = _DirichletData(_mesh.coordinates(index1-1,0), _mesh.coordinates(index1-1,1));
         	
       		//Set Dirichlet data in solution vector:
       		fl_u(index1) = fl_uD(index1);
       		
        }
    }
    
    for (int i=0; i<_mesh.coupling.crossPointsBdryData.length(); ++i) {
    
        if (_mesh.coupling.crossPointsBdryData(i)!=0) {
        
            //Set Dirichlet data on nodes:
      		fl_uD(i+1) = _DirichletData(_mesh.coordinates(i,0), _mesh.coordinates(i,1));
      		
      		//Set Dirichlet data in solution vector:
      		fl_u(i+1) = fl_uD(i+1);
      		
        }
    }  
      
}    


