#ifndef GS_MPI_BLAS_CPP
#define GS_MPI_BLAS_CPP 1

#include "gs_mpi_blas.h"


//FLENS-based dense GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max())
{
	using namespace flens;

    typedef typename VB::ElementType  		ElementType;
    typedef typename VB::IndexType    		IndexType;
    typedef VB				    			VectorType;
    typedef DenseVector<Array<IndexType> >	IVector;

    const Underscore<IndexType> _; // FLENS operator (range access)

    const Coupling &coupling = b.coupling;

  	// GS_MPI starts here
  	// First resort indices
    IndexType numNodes=x.length();
    
    /* Step 0.1:
     *    Build index vectors for 
     *      - cross points (indexV), 
     *      - boundary nodes (indexE),
     *      - inner nodes (indexI) 
     */
    
    // cross point are listed first
    int nV=coupling.local2globalCrossPoints.length();
    IVector indexV(nV);
    for(IndexType k=1; k<=nV; ++k)
    {
        indexV(k) = k;
    }
    
    // build innerNodes with innerNodes(k)==0 => k-th node  is an inner node 
    IVector innerNodes(numNodes);
    for(IndexType k=1; k<=nV; ++k)
    {
         innerNodes(k) = 1;
    }
    IndexType nE=0, counter=1;
    
    // count number of boundary nodes
    for(int k=0; k<coupling.numCoupling; ++k)
    {
         nE += coupling.boundaryNodes[k].length()-2;
    }
    IVector indexE(nE);
    for(int k=0; k<coupling.numCoupling; ++k)
    { 
        for(int j=2; j<=coupling.boundaryNodes[k].length()-1; ++j)
        {
            indexE(counter) =  coupling.boundaryNodes[k](j);
            innerNodes(coupling.boundaryNodes[k](j)) = 1;
            counter ++;
        }
    }
    
    // set innerNodes for fixed Nodes (Dirichlet nodes)
    for(int k=1; k<=bc.length(); ++k)
    {
        innerNodes(bc(k)) = 1;
    }
    // count inner nodes
    IndexType nI=0;
    for(IndexType j=nV+1; j<=numNodes; ++j)
    {
        if(innerNodes(j)==0) ++nI;
    }

        
    IVector indexI(nI);
    counter=1;
    for(IndexType j=nV+1; j<=numNodes; ++j)
    {
        if(innerNodes(j)==0)
        {
            indexI(counter) = j;
            counter++;
        }
    }
    
    /* Step 0.2:
     *    Set diagonal block matrices A_VV and A_EE 
     *    of matrix A
     *      - first local values
     *      - communication to get global values
     */
    
    // set local values for A_EE and A_VV
    FLENSDataVector diag(numNodes, coupling, flens::typeII);
    VectorType subdiag(nE);
    IndexType idxk;
    // set values at diagonal
    for(IndexType k=1; k<=numNodes; ++k)
    {
        diag(k) = A(k,k);
    }

    // Now global values: 
    // a) communication for diagonal 
    diag.typeII_2_I();
    // b) communication for subdiagonal
    int offset=0;
    for (int j=0; j<coupling.numCoupling; ++j)
    {
        for (int i=1; i<=coupling.maxColor; ++i)
        {
            if (coupling.colors(j+1) == i){
                IndexType numBdryNodes = coupling.boundaryNodes[j].length()-2;
                // only communicate if there are more than 1 boundary nodes on coupling boundary (no cross Points!)
                if(numBdryNodes>1){
  
                    IVector sendIndex = coupling.boundaryNodes[j](_(2,numBdryNodes+1));
                  
                    VectorType u_send1(numBdryNodes-1);
                    VectorType u_recv1(numBdryNodes-1);
                    // set local values
                    for(IndexType k=1; k<=numBdryNodes-1; ++k)
                    {
                        u_send1(k) =  A(sendIndex(k+1), sendIndex(k));
                    }
                    // get values from other processes
                    MPI::COMM_WORLD.Sendrecv(u_send1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j+1)-1, 0,
                                   u_recv1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j+1)-1, 0);
                    // add values from other processes (!! numbering is opposite !!)
                    for(IndexType k=1; k<=numBdryNodes-1; ++k)
                    {
                        subdiag(k+offset) =  u_send1(k) + u_recv1(numBdryNodes-k);
                    }
                }
                offset+=numBdryNodes;
            }
        }
    }    
    // extract matrices A_VV and  A_EE
    VectorType A_VV(nV), A_EE_diag(nE), A_EE_udiag(nE), A_EE_ldiag(nE);
    for(IndexType k=1; k<=nV; ++k)
    {
        A_VV(k) = diag(k);
    }
    for(IndexType k=1; k<=nE; ++k)
    {
        idxk=indexE(k);
        A_EE_diag(k) = diag(idxk);
    }
    /* 
     * Start iteration  
     */
     
    // initialize residual 
    FLENSDataVector r(numNodes, coupling, flens::nonMPI);
    // set x to zero at fixed nodes
    for(int i=1; i<=bc.length(); ++i) {
          x(bc(i)) = 0;
    }

    // set bdryNodes vector as we only solve at free nodes! 
    IVector bdryNodes(numNodes-nI);
    for(IndexType k=1; k<=nV; ++k) bdryNodes(k) = k;
    counter=nV+1;
    for(IndexType k=nV+1; k<=numNodes; ++k)
    {
        if(innerNodes(k)!=0)
        { 
            bdryNodes(counter) = k;
            ++counter;
        }
    }

    for(int ell=1; ell<=maxIterations; ++ell)
    {
        
        /* Step 1.1:
         *    Calculate residual on crosspoints
         *    r_V = b_V - A_VV*x_V - A_VE*x_E - A_VI * x_I
         */
        for(IndexType k=1; k<=nV; ++k)
        {
            if(coupling.crossPointsBdryData(k)==0)
            {
                ElementType tmp2=0;
                for(IndexType j=1; j<=A.lastCol(); ++j)
                { 
                    tmp2 += A(k,j)*x(j);
                }    
                r(k) =  b(k) - tmp2;
            }        
        }
        /* Step 1.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_V),
         *    saved in r as well
         */
        r.commCrossPoints();
        
        /* Step 1.3:
         *    Update x at cross points 
         */
        for(IndexType k=1; k<=nV; ++k)
        {
            x(k) += r(k)/diag(k);
        }  
        /* Step 2.1:
         *    Calculate residual on boundary nodes 
         *    r_E = b_E - A_EV*x_V - A_EE*x_E - A_EI * x_I
         */
        for(IndexType k=1; k<=nE; ++k)
        {
            IndexType idx = indexE(k);
            ElementType tmp2=0;
            for(IndexType j=1; j<=A.lastCol(); ++j)
            {
                tmp2 += A(idx,j)*x(j);
            }    
            r(idx) = b(idx) - tmp2;    
        }
        
        /* Step 2.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_E),
         *    saved in r as well
         */
        r.commBoundaryNodes();
        
        /* Step 2.3:
         *    Update x at boundary nodes:
         *      - compute A_EE^-1*rE
         *      - update x
         */
          
        // Computation of A_EE^-1*rE:
        VectorType xE(nE), rE(nE);
        //  a) copy subdiagonals (as they are changed in solveTridiag)

		A_EE_ldiag(_(1,nE-1)) = subdiag(_(1,nE-1));
        A_EE_udiag(_(1,nE-1)) = subdiag(_(1,nE-1));

        //  b) set right hand side
        for(IndexType k=1; k<=nE; ++k)
        {
             rE(k) =  r(indexE(k));
        }
        //  c) invoke tridiagonal solver
        solveTridiag(A_EE_ldiag, A_EE_diag, A_EE_udiag,xE,rE);
        
        // Update x at boundary nodes
        for(IndexType k=1; k<=nE; ++k)
        { 
            x(indexE(k)) += xE(k);    
        }
        /* Step 3:
         *    Calculate x on inner nodes: 
         *    Gauß-Seidel step in forward direction (only at FREE NODES)
         */
        for(IndexType j=1; j<=numNodes; ++j)
        {
            if(innerNodes(j)==0)
            {    
                ElementType tmp1=0;
                for(IndexType mu=1; mu<=A.lastCol(); ++mu)
                {
                    tmp1 += A(j,mu)*x(mu);
                }
        
                x(j) += (b(j)-tmp1)/A(j,j);
            }
        }
    }
    
    return maxIterations;
};


//FLENS-based sparse GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_mpi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max())
{
    using namespace flens;

    typedef typename VB::ElementType                ElementType;
    typedef typename VB::IndexType                  IndexType;
    typedef VB                                      VectorType;
    typedef DenseVector<Array<IndexType> >          IVector;

    const Underscore<IndexType> _; // FLENS operator (range access)

    ElementType Zero(0);
    const Coupling &coupling = b.coupling;

    // GS_MPI starts here
    // First resort indices
    IndexType numNodes=x.length();
    
    /* Step 0.1:
     *    Build index vectors for 
     *      - cross points (indexV), 
     *      - boundary nodes (indexE),
     *      - inner nodes (indexI) 
     */
    
    // cross point are listed first
    int nV=coupling.local2globalCrossPoints.length();
    IVector indexV(nV);
    for(IndexType k=1; k<=nV; ++k)
    {
        indexV(k) = k;
    }
    
    // build innerNodes with innerNodes(k)==0 => k-th node  is an inner node 
    IVector innerNodes(numNodes);
    for(IndexType k=1; k<=nV; ++k)
    {
         innerNodes(k) = 1;
    }
    IndexType nE=0, counter=1;
    
    // count number of boundary nodes
    for(int k=0; k<coupling.numCoupling; ++k)
    {
         nE += coupling.boundaryNodes[k].length()-2;
    }
    IVector indexE(nE);
    for(int k=0; k<coupling.numCoupling; ++k)
    { 
        for(int j=2; j<=coupling.boundaryNodes[k].length()-1; ++j)
        {
            indexE(counter) =  coupling.boundaryNodes[k](j);
            innerNodes(coupling.boundaryNodes[k](j)) = 1;
            counter ++;
        }
    }

    // set innerNodes for fixed Nodes (Dirichlet nodes)
    for(int k=1; k<=bc.length(); ++k)
    {
        innerNodes(bc(k)) = 1;
    }
    // count inner nodes
    IndexType nI=0;
    for(IndexType j=nV+1; j<=numNodes; ++j)
    {
        if(innerNodes(j)==0) ++nI;
    }

        
    IVector indexI(nI);
    counter=1;
    for(IndexType j=nV+1; j<=numNodes; ++j)
    {
        if(innerNodes(j)==0)
        {
            indexI(counter) = j;
            counter++;
        }
    }
    
    /* Step 0.2:
     *    Set diagonal block matrices A_VV and A_EE 
     *    of matrix A
     *      - first local values
     *      - communication to get global values
     */
    
    // set local values for A_EE and A_VV
    FLENSDataVector diag(numNodes, coupling, flens::typeII);
    VectorType subdiag(nE);
    IndexType idxk;

    // Access values of A
    const auto &_rows = A.engine().rows();
    const auto &_cols = A.engine().cols();
    const auto &_vals = A.engine().values();

    ElementType Akk = Zero;
    // set values at diagonal
    for(IndexType k=1; k<=numNodes; ++k)
    {
        // Locate value
        for (IndexType j=_rows(k); j<_rows(k+1); ++j)
        {
            if (_cols(j)==k)
            {
                Akk = _vals(j);
            }
        }
        assert(Akk!=Zero);
        diag(k) = Akk;
    }


    // Now global values: 
    // a) communication for diagonal 
    diag.typeII_2_I();
    // b) communication for subdiagonal
    int offset=0;
    for (int j=0; j<coupling.numCoupling; ++j)
    {
        for (int i=1; i<=coupling.maxColor; ++i)
        {
            if (coupling.colors(j+1) == i){
                IndexType numBdryNodes = coupling.boundaryNodes[j].length()-2;
                // only communicate if there are more than 1 boundary nodes on coupling boundary (no cross Points!)
                if(numBdryNodes>1){
                    
                    IVector sendIndex = coupling.boundaryNodes[j](_(2,numBdryNodes+1));
                    
                    VectorType u_send1(numBdryNodes-1);
                    VectorType u_recv1(numBdryNodes-1);
                    // set local values
                    for(IndexType k=1; k<=numBdryNodes-1; ++k)
                    {
                        Akk = Zero;
                        // Locate value
                        for (IndexType it=_rows(sendIndex(k+1)); it<_rows(sendIndex(k+1)+1); ++it)
                        {
                            if (_cols(it)==sendIndex(k))
                            {
                                Akk = _vals(it);
                            }
                        }
                        u_send1(k) =  Akk;
                    }
                    // get values from other processes
                    MPI::COMM_WORLD.Sendrecv(u_send1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j+1)-1, 0,
                                   u_recv1.data() , numBdryNodes-1 , MPI::DOUBLE,
                                   coupling.neighbourProcs(j+1)-1, 0);
                    // add values from other processes (!! numbering is opposite !!)
                    for(IndexType k=1; k<=numBdryNodes-1; ++k)
                    {
                        subdiag(k+offset) =  u_send1(k) + u_recv1(numBdryNodes-k);
                    }
                }
                offset+=numBdryNodes;
            }
        }
    }    
    // extract matrices A_VV and  A_EE
    VectorType A_VV(nV), A_EE_diag(nE), A_EE_udiag(nE), A_EE_ldiag(nE);
    for(IndexType k=1; k<=nV; ++k)
    {
        A_VV(k) = diag(k);
    }
    for(IndexType k=1; k<=nE; ++k)
    {
        idxk=indexE(k);
        A_EE_diag(k) = diag(idxk);
    }
    /* 
     * Start iteration  
     */
     
    // initialize residual 
    FLENSDataVector r(numNodes, coupling, flens::nonMPI);
    // set x to zero at fixed nodes
    for(int i=1; i<=bc.length(); ++i) {
          x(bc(i)) = 0;
    }

    // set bdryNodes vector as we only solve at free nodes! 
    IVector bdryNodes(numNodes-nI);
    for(IndexType k=1; k<=nV; ++k) bdryNodes(k) = k;
    counter=nV+1;
    for(IndexType k=nV+1; k<=numNodes; ++k)
    {
        if(innerNodes(k)!=0)
        { 
            bdryNodes(counter) = k;
            ++counter;
        }
    }

    for(int ell=1; ell<=maxIterations; ++ell)
    {
        
        /* Step 1.1:
         *    Calculate residual on crosspoints
         *    r_V = b_V - A_VV*x_V - A_VE*x_E - A_VI * x_I
         */
        for(IndexType k=1; k<=nV; ++k)
        {
            if(coupling.crossPointsBdryData(k)==0)
            {
                ElementType tmp2=0;
                for(IndexType j=_rows(k); j<_rows(k+1); ++j)
                { 
                    tmp2 += _vals(j)*x(_cols(j));
                }    
                r(k) =  b(k) - tmp2;
            }        
        }
        /* Step 1.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_V),
         *    saved in r as well
         */
        r.commCrossPoints();
        
        /* Step 1.3:
         *    Update x at cross points 
         */
        for(IndexType k=1; k<=nV; ++k)
        {
            x(k) += r(k)/diag(k);
        }  
        /* Step 2.1:
         *    Calculate residual on boundary nodes 
         *    r_E = b_E - A_EV*x_V - A_EE*x_E - A_EI * x_I
         */
        for(IndexType k=1; k<=nE; ++k)
        {
            IndexType idx = indexE(k);
            ElementType tmp2=0;
            for(IndexType j=_rows(idx); j<_rows(idx+1); ++j)
            {
                tmp2 += _vals(j)*x(_cols(j));
            }    
            r(idx) = b(idx) - tmp2;    
        }
        
        /* Step 2.2:
         *    Communicate with other proceses
         *    to compute typeI residual (w_E),
         *    saved in r as well
         */
        r.commBoundaryNodes();
        
        /* Step 2.3:
         *    Update x at boundary nodes:
         *      - compute A_EE^-1*rE
         *      - update x
         */
          
        // Computation of A_EE^-1*rE:
        VectorType xE(nE), rE(nE);
        //  a) copy subdiagonals (as they are changed in solveTridiag)

        A_EE_ldiag(_(1,nE-1)) = subdiag(_(1,nE-1));
        A_EE_udiag(_(1,nE-1)) = subdiag(_(1,nE-1));

        //  b) set right hand side
        for(IndexType k=1; k<=nE; ++k)
        {
             rE(k) =  r(indexE(k));
        }
        //  c) invoke tridiagonal solver
        solveTridiag(A_EE_ldiag, A_EE_diag, A_EE_udiag,xE,rE);
        
        // Update x at boundary nodes
        for(IndexType k=1; k<=nE; ++k)
        { 
            x(indexE(k)) += xE(k);    
        }
        /* Step 3:
         *    Calculate x on inner nodes: 
         *    Gauß-Seidel step in forward direction (only at FREE NODES)
         */
        for(IndexType j=1; j<=numNodes; ++j)
        {
            if(innerNodes(j)==0)
            {    
                ElementType tmp1=0;
                for(IndexType mu=_rows(j); mu<_rows(j+1); ++mu)
                {
                    tmp1 += _vals(mu)*x(_cols(mu));
                }
        
                x(j) += (b(j)-tmp1)/diag(j);
            }
        }
    }
    
    return maxIterations;
};
   

template <typename V>
void
solveTridiag(V &ldiag, V &diag, V &udiag, V &x, V &b)
{
	using namespace flens;
	typedef typename V::ElementType  		ElementType;
    typedef typename V::IndexType    		IndexType;

	/* *** Thomas algorithm to solve tridiagonal system */ 
	IndexType n=diag.length();
	if(n==1) x(1) = b(1)/diag(1);
	else
	{
		/* *** forward iteration */
		udiag(1) = udiag(1)/diag(1);
		b(1) = b(1)/diag(1);
		for(IndexType k=2; k<=n-1; ++k)
		{
			udiag(k) = udiag(k)/( diag(k)-udiag(k-1)*ldiag(k-1) );
			b(k) = ( b(k)-b(k-1)*ldiag(k-1) )/( diag(k)-udiag(k-1)*ldiag(k-1) );
		}
		b(n) = (b(n)-b(n-1)*ldiag(n-1))/(diag(n)-udiag(n-1)*ldiag(n-1));
		/* *** backward iteration */
		x(n) = b(n);
		for(int k=n-1; k>0; --k) x(k) = b(k)-udiag(k)*x(k+1);
	}
};

#endif	//GS_MPI_BLAS_CPP