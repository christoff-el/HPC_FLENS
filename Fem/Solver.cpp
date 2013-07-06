#include "Solver.hpp"

/********************************************************************************************************************/
/**********************************************     CG - Method      ************************************************/
/********************************************************************************************************************/

int CG(CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes,  int maxIt, double tol)
{
    
    //Solve for x using CG method under FLENS framework:
    maxIt = cg_nompi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);

    return maxIt;    
}

int CG_MPI(CRSMatrix &A, DataVector &x, DataVector &b, IndexVector &dirichletNodes,  int maxIt, double tol)
{
	
	//Solve for x using MPI CG method under FLENS framework:
	maxIt = cg_mpi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);
	
  	return maxIt;    
}

/********************************************************************************************************************/
/*******************************************   Gauß-Seidel - Method    **********************************************/
/********************************************************************************************************************/

int forwardGS( CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes, int maxIt, double tol)
{
    maxIt = gs_dense_nompi_blas_wrapper(A, x, b, dirichletNodes, maxIt, tol);
        
    return maxIt;    
}

int forwardGS_MPI( CRSMatrix &A, DataVector &x, DataVector &b, Coupling &coupling, 
                   IndexVector &dirichletNodes, int maxIt)
{

    assert(A.numRows()==b.values.length() && A.numCols()==x.values.length());
    assert(x.type==typeI && b.type==typeII);
    
    int numNodes=x.values.length();
    
    /* Step 0.1:
     *    Build index vectors for 
     *      - cross points (indexV), 
     *      - boundary nodes (indexE),
     *      - inner nodes (indexI) 
     */
    
    // cross point are listed first
    int nV=coupling.local2globalCrossPoints.length();
    IndexVector indexV(nV);
    for(int k=0;k<nV;k++){
        indexV(k) = k+1;
    }
    
    // build innerNodes with innerNodes(k)==0 => k-th node  is an inner node 
    IndexVector innerNodes(numNodes);
    for(int k=0;k<nV;k++){
         innerNodes(k) = 1;
    }
    int nE=0, counter=0;
    
    // count number of boundary nodes
    for(int k=0;k<coupling.numCoupling;k++){
         nE += coupling.boundaryNodes[k].length()-2;
    }
    IndexVector indexE(nE);
    for(int k=0;k<coupling.numCoupling;k++){ 
        for(int j=1;j<coupling.boundaryNodes[k].length()-1;j++){
            indexE(counter) =  coupling.boundaryNodes[k](j);
            innerNodes(coupling.boundaryNodes[k](j)-1) = 1;
            counter ++;
        }
    }
    
    // set innerNodes for fixed Nodes (Dirichlet nodes)
    for(int k=0;k<dirichletNodes.length();k++) {
        innerNodes(dirichletNodes(k)-1) = 1;
    }
    // count inner nodes
    int nI=0;
    for(int j=nV;j<numNodes;j++){
        if(innerNodes(j)==0) nI++;
    }
        
    IndexVector indexI(nI);
    counter=0;
    for(int j=nV;j<numNodes; j++){
        if(innerNodes(j)==0){
            indexI(counter) = j+1;
            counter++;
        }
    }
    
    /* Step 0.2:
     *    Set diagonal block matrices A_VV and A_EE 
     *    of matrix A
     *      - first local values
     *      - communication to get global values
     */
    double * Adata = A.data();
    int * ArowPtr = A.rowPtr();
    int * AcolIndex = A.colIndex();
    
    // set local values for A_EE and A_VV
    DataVector diag(coupling,numNodes,typeII);
    Vector subdiag(nE);
    int idxk;
    // set values at diagonal
    for(int k=0;k<numNodes;k++){
        diag.values(k) = A(k,k);
    }

    // Now global values: 
    // a) communication for diagonal 
    diag.typeII_2_typeI();
    // b) communication for subdiagonal
    int offset=0;
    for (int j=0; j<coupling.numCoupling; j++) {
        for (int i=1; i<=coupling.maxColor; i++) {
            if (coupling.colors(j) == i){
                int numBdryNodes = coupling.boundaryNodes[j].length()-2;
                // only communicate if there are more than 1 boundary nodes on coupling boundary (no cross Points!)
                if(numBdryNodes >1){
                    IndexVector sendIndex(numBdryNodes);
                    sendIndex.set(0,numBdryNodes, coupling.boundaryNodes[j].data()+1);
                    Vector u_send1(numBdryNodes-1);
                    Vector u_recv1(numBdryNodes-1);
                    // set local values
                    for( int k=0; k<numBdryNodes-1; k++) {
                        u_send1(k) =  A(sendIndex(k+1)-1, sendIndex(k)-1);
                    }
                    // get values from other processes
                    MPI::COMM_WORLD.Sendrecv(u_send1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j)-1, 0,
                                   u_recv1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j)-1, 0);
                    // add values from other processes (!! numbering is opposite !!)
                    for( int k=0; k<numBdryNodes-1; k++) {
                        subdiag(k+offset) =  u_send1(k) + u_recv1(numBdryNodes-k-2);
                    }
                }
                offset+=numBdryNodes;
            }
        }
    }    
    // extract matrices A_VV and  A_EE
    Vector A_VV(nV), A_EE_diag(nE), A_EE_udiag(nE),A_EE_ldiag(nE );
    for(int k=0;k<nV;k++){
        A_VV(k) = diag.values(k);
    }
    for(int k=0; k<nE;k++){
        idxk=indexE(k)-1;
        A_EE_diag(k) = diag.values(idxk);
    }
    
    /* 
     * Start iteration  
     */
     
    // initialize residual 
    DataVector r(coupling,numNodes);
    // set x to zero at fixed nodes
    for(int i = 0; i < dirichletNodes.length(); i++) {
          x.values( dirichletNodes(i)-1 ) = 0;
    }
    // set bdryNodes vector as we only solve at free nodes! 
    IndexVector bdryNodes(numNodes-nI);
    for(int k=0;k<nV;k++) bdryNodes(k) = k+1;
    counter=nV;
    for(int k=nV; k<numNodes; k++){
        if(innerNodes(k)!=0){ 
            bdryNodes(counter) = k+1;
            counter ++;
        }
    }
    for(int ell=0; ell<maxIt; ell++){
        
        /* Step 1.1:
         *    Calculate residual on crosspoints
         *    r_V = b_V - A_VV*x_V - A_VE*x_E - A_VI * x_I
         */
        for(int k=0; k< nV; k++){
            if(coupling.crossPointsBdryData(k)==0){
                double tmp2=0;
                for(int j=ArowPtr[k];j<ArowPtr[k+1];j++){ 
                    tmp2 += Adata[j]*x.values(AcolIndex[j]);
                }    
                r.values(k) =  b.values(k) - tmp2;
            }        
        }
        /* Step 1.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_V),
         *    saved in r as well
         */
        r.communicationCrossPoints();
        
        /* Step 1.3:
         *    Update x at cross points 
         */
        for(int k=0;k<nV;k++){
            x.values(k) += r.values(k)/A_VV(k);
        }
        
        /* Step 2.1:
         *    Calculate residual on boundary nodes 
         *    r_E = b_E - A_EV*x_V - A_EE*x_E - A_EI * x_I
         */
        for(int k=0; k<nE; k++){
            int idx = indexE(k)-1;
            double tmp2=0;
            for(int j=ArowPtr[idx];j<ArowPtr[idx+1];j++){
                tmp2 += Adata[j]*x.values(AcolIndex[j]);
            }    
            r.values(idx) = b.values(idx) - tmp2;    
        }
        
        /* Step 2.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_E),
         *    saved in r as well
         */
        r.communicationBoundaryNodes();
        
        /* Step 2.3:
         *    Update x at boundary nodes:
         *      - compute A_EE^-1*rE
         *      - update x
         */
          
        // Computation of A_EE^-1*rE:
        Vector xE(nE), rE(nE);
        //  a) copy subdiagonals (as they are changed in solveTridiag)
        A_EE_ldiag.set(0,nE-1,subdiag.data());
        A_EE_udiag.set(0,nE-1,subdiag.data());
        //  b) set right hand side
        for(int k=0; k<nE; k++){
             rE(k) =  r.values(indexE(k)-1);
        }
        //  c) invoke tridiagonal solver
        solveTridiag(A_EE_ldiag, A_EE_diag, A_EE_udiag,xE,rE);
        
        // Update x at boundary nodes
        for(int k=0;k<nE;k++){ 
            x.values(indexE(k)-1) += xE(k);    
        }    
        
        /* Step 3:
         *    Calculate x on inner nodes: 
         *    Gauß-Seidel step in forward direction (only at FREE NODES)
         */
        for(int j=0;j<numNodes; j++){
            if(innerNodes(j)==0){    
                double tmp1=0;
                for(int mu=ArowPtr[j];mu<ArowPtr[j+1];mu++){
                    tmp1 += Adata[mu]*x.values(AcolIndex[mu]);
                }
        
                x.values(j) += ( b.values(j)-tmp1 )/A(j,j);
            }
        }
    }
    return maxIt;
}//end GS


void solveTridiag(Vector &ldiag, Vector &diag, Vector &udiag, Vector &x, Vector &b)
{
	/* *** Thomas algorithm to solve tridiagonal system */ 
	int n=diag.length();
	if(n==1) x(0) = b(0)/diag(0);
	else{
		/* *** forward iteration */
		udiag(0) = udiag(0)/diag(0);
		b(0) = b(0)/diag(0);
		for(int k=1;k<n-1;k++){
			udiag(k) = udiag(k)/( diag(k)-udiag(k-1)*ldiag(k-1) );
			b(k) = ( b(k)-b(k-1)*ldiag(k-1) )/( diag(k)-udiag(k-1)*ldiag(k-1) );
		}
		b(n-1) = ( b(n-1)-b(n-2)*ldiag(n-2) )/( diag(n-1)-udiag(n-2)*ldiag(n-2) );
		/* *** backward iteration */
		x(n-1) = b(n-1);
		for(int k=n-2;k>=0;k--) x(k) = b(k)-udiag(k)*x(k+1);
	}
}
