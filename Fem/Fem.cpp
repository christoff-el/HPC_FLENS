#ifndef FEM_CPP
#define FEM_CPP 1

#include "Fem.hpp"

/************************************************************************************************************/
/************************************    constructor   ******************************************************/
/************************************************************************************************************/

//'Default' constructor when an unsupported method specialisation instantiated:
template <typename METH>
FEM<METH>::FEM(Mesh &mesh, double (*f)(double,double), 
        			double (*DirichletData)(double,double), double (*g)(double,double))
{
	METH::CHK;
}

//Non-MPI constructor:
template <>
FEM<flens::MethNonMPI>::FEM(Mesh &mesh, double (*f)(double,double), 
        			double (*DirichletData)(double,double), double (*g)(double,double))
        :	_mesh(mesh),
        	_A(),
	    	_uD(mesh.numNodes),
         	_u(mesh.numNodes),
         	_b(mesh.numNodes),
         	_f(f),
         	_g(g),
         	_DirichletData(DirichletData)
{
	//Constructor to initialise FLENSDataVector<FLvNonMPI>'s without specifying coupling.
}

//MPI constructor:
template <>
FEM<flens::MethMPI>::FEM(Mesh &mesh, double (*f)(double,double), 
        			double (*DirichletData)(double,double), double (*g)(double,double))
        :	_mesh(mesh),
        	_A(),
	    	_uD(mesh.numNodes, mesh.coupling),
         	_u(mesh.numNodes, mesh.coupling),
         	_b(mesh.numNodes, mesh.coupling),
         	_f(f),
         	_g(g),
         	_DirichletData(DirichletData)
{
	//Constructor to initialise FLENSDataVector<FLvTypeI/II>s> specifying coupling.
}

template <typename METH>
FEM<METH>::~FEM()
{
    
}


/**********************************************************************************************/
/*******************************    mehtods for assembling   **********************************/
/**********************************************************************************************/
// assemble the Galerkin matrix and the right hand side 
template <typename METH>
void 
FEM<METH>::assemble()
{
  	//Write values of dirichlet data:
    _updateDirichlet();
    
    flens::DenseVector<flens::Array<int> > 	I(9*_mesh.numElements);
    flens::DenseVector<flens::Array<int> > 	J(9*_mesh.numElements);
    flens::DenseVector<flens::Array<double> > val(9*_mesh.numElements);
    double c1[2], d21[2], d31[2];
    double area4, fval, a, b, c;
    
    for (int i=0; i<_mesh.numElements; ++i)
    {
  
        /*** Assemble stiffness matrix A ***/
        
        //Compute area of element: 
        // c1 = first vertex of element, 
        // d21 vector from first to second vertex, 
        // d31 vector from third to first vertex:
        c1[0]  = _mesh.coordinates( _mesh.elements(i+1,1),1 ); //<-----numbering starts at 1
        c1[1]  = _mesh.coordinates( _mesh.elements(i+1,1),2 );
        d21[0] = _mesh.coordinates( _mesh.elements(i+1,2),1 ) - c1[0];
        d21[1] = _mesh.coordinates( _mesh.elements(i+1,2),2 ) - c1[1];
        d31[0] = _mesh.coordinates( _mesh.elements(i+1,3),1 ) - c1[0];
        d31[1] = _mesh.coordinates( _mesh.elements(i+1,3),2 ) - c1[1];
        
        //Compute 4*|T_i| via 2-dimensional cross product:
        area4 = 2*(d21[0]*d31[1] - d21[1]*d31[0]);
        
        //Set I and J  (row and column index)
        for(int j=0; j<3; ++j) 
        {
        
            I(i*9+j  +1) = _mesh.elements(i+1,j+1);
            I(i*9+j+3+1) = _mesh.elements(i+1,j+1);
            I(i*9+j+6+1) = _mesh.elements(i+1,j+1);
        
            J(i*9+3*j  +1) = _mesh.elements(i+1,j+1);
            J(i*9+3*j+1+1) = _mesh.elements(i+1,j+1);
            J(i*9+3*j+2+1) = _mesh.elements(i+1,j+1);
            
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
        _b( _mesh.elements(i+1,1)) += fval;
        _b( _mesh.elements(i+1,2)) += fval;
        _b( _mesh.elements(i+1,3)) += fval;
        
    }

    
    /*** Incorprate Neumann data ***/

    //Compute bj += int _{Gamma_N} g(x) phi_j(x) dx:
    for(int k=0; k<_mesh.numNeumann; ++k) {
    
        //Walk through all neumann edges:
        for(int j=1; j<=_mesh.neumann[k].length()-1; ++j) {
        
            double cn1[2],cn2[2], length2;
            cn1[0]  = _mesh.coordinates(_mesh.neumann[k](j),1);
            cn1[1]  = _mesh.coordinates(_mesh.neumann[k](j),2);
            cn2[0]  = _mesh.coordinates(_mesh.neumann[k](j+1),1);
            cn2[1]  = _mesh.coordinates(_mesh.neumann[k](j+1),2);
            length2 = sqrt( (cn1[0]-cn2[0])*(cn1[0]-cn2[0]) + (cn1[1]-cn2[1])*(cn1[1]-cn2[1]) ) / 2.;
            
            //Evaluate Neumann data g at midpoint of neumann edge and multiply with half length:
            double gmE = length2*_g(0.5*(cn1[0]+cn2[0]) , 0.5*(cn1[1]+cn2[1]));
            _b( _mesh.neumann[k](j)  ) += gmE;
            _b( _mesh.neumann[k](j+1)) += gmE;   
             
        }   
         
    }


    /*** Set stiffness matrix ***/
    
    //Build FLENS CRS matrix from I, J, vals (rows, cols, values):
    flens::GeCoordMatrix<flens::CoordStorage<double> > fl_A_coord(_uD.length(),_uD.length());

    for (int i=1; i<=I.length(); ++i) {

        if (val(i)!=0) {
            fl_A_coord(I(i),J(i)) += val(i);
        }
    	
    }

    _A = fl_A_coord;


    /*** Set right-hand side vector b ***/    
    flens::FLENSDataVector<flens::FLvTypeII> fl_Au(_uD.length(), _mesh.coupling);
        

    /*** Incorporate Dirichlet-Data ***/
    flens::blas::mv(flens::NoTrans, 1., _A, _uD, 0., fl_Au);
    flens::blas::axpy(-1., fl_Au, _b);    

}


/***********************************************************************************************************/
/**********************************    methods for solving SLE   *******************************************/
/***********************************************************************************************************/
template <typename METH>
void
FEM<METH>::solve(Solver method)
{    

    /*** Set fixed nodes (all Dirichletnodes) ***/
    
    //Count the number of fixed nodes:
    int numfixed = 0;
    for (int k=0; k<_mesh.numDirichlet; ++k) {
        numfixed += _mesh.dirichlet[k].length();
    } 
    
    for (int k=1; k<=_mesh.coupling.crossPointsBdryData.length(); ++k) {
        if(_mesh.coupling.crossPointsBdryData(k)!=0) numfixed++;        
    }
    
    //Set fixedNodes:
    flens::DenseVector<flens::Array<int> > fixedNodes(numfixed);
    int index=1;
    for (int k=0; k<_mesh.numDirichlet; ++k) { 
    
    	//Setting:
    	for (int l=index; l<index+_mesh.dirichlet[k].length(); ++l) {
    		fixedNodes(l) = _mesh.dirichlet[k](l);			//fixedNodes.set(index, _mesh.dirichlet[k].length(), _mesh.dirichlet[k].data() );
        }
        
        index += _mesh.dirichlet[k].length();
        
    }
    
    for (int k=1; k<=_mesh.coupling.crossPointsBdryData.length(); ++k) {
    
        if (_mesh.coupling.crossPointsBdryData(k) != 0) {
        
            fixedNodes(index) = k;
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
        //if (!_mesh.distributed()) {
        if (std::is_same<METH, flens::MethNonMPI>::value) {
            it = cg_nompi_blas(_A ,_b, _u, fixedNodes, maxIt, tol);
        }
        
        //Parallel solver:
        else { 
            it = cg_mpi_blas(_A ,_b, _u, fixedNodes,  maxIt, tol);
        }
        
    }
    
    //GS Solver:
    else if (method == gs) {
    	
    	
    	//Serial solver:
		if (std::is_same<METH, flens::MethNonMPI>::value) {
			it = gs_nompi_blas(_A, _b, _u, fixedNodes, maxIt, tol);
		}
		
		//Parallel solver:
		else {
			it = gs_mpi_blas(_A, _b, _u, fixedNodes, maxIt);
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
template <typename METH>
void 
FEM<METH>::refineRed()
{
  /* In future projects this function will be modified such that all meshes
     are stored in a hierarchy object which is needed for the multigrid algorithm.
     For now we don't need this and the function just refines the mesh and
     adjusts the size of the vectors 
  */
  
    /*** Refine mesh ***/
	_mesh.refineRed();

    /*** Resize DataVectors (note: clears data) ***/
  	_u.resize(_mesh.numNodes);
  	_b.resize(_mesh.numNodes);
  	_uD.resize(_mesh.numNodes);
        
}

/**********************************************************************************************/
/*******************************  getter and write methods   **********************************/
/**********************************************************************************************/
template <typename METH>
flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1> > >
FEM<METH>::getA()
{

  	return _A;

}

template <typename METH>
typename SelectDataVector<METH>::TypeII
FEM<METH>::getb()
{
	
    return _b;
    
}

template <typename METH>
int 
FEM<METH>::getNumElements()
{
    return _mesh.numElements;
}

template <typename METH>
void 
FEM<METH>::writeSolution(int proc, std::string filename)
{
	_mesh.writeData(proc, filename);
	_u.writeData(proc, filename + "solution");
}

/**********************************************************************************************/
/*********************************    private methods   ***************************************/
/**********************************************************************************************/

/* write the dirichlet data into the vector _uD and in the solution _u */
template <typename METH>
void 
FEM<METH>::_updateDirichlet()
{

    int index1;
    for (int k=0; k<_mesh.numDirichlet; ++k) {
        for (int i=1; i<=_mesh.dirichlet[k].length(); ++i) {
        
        	index1 = _mesh.dirichlet[k](i);
        	
            //Set Dirichlet data on nodes:
         	_uD(index1) = _DirichletData(_mesh.coordinates(index1,1), _mesh.coordinates(index1,2));
         	
       		//Set Dirichlet data in solution vector:
       		_u(index1) = _uD(index1);
       		
        }
    }
    
    for (int i=1; i<=_mesh.coupling.crossPointsBdryData.length(); ++i) {
    
        if (_mesh.coupling.crossPointsBdryData(i)!=0) {
        
            //Set Dirichlet data on nodes:
      		_uD(i) = _DirichletData(_mesh.coordinates(i,1), _mesh.coordinates(i,2));
      		
      		//Set Dirichlet data in solution vector:
      		_u(i) = _uD(i);
      		
        }
    }  
      
}    


#endif	//FEM_CPP
