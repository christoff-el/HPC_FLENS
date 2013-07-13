#include "Mesh.hpp"

// for RGB-refinement
enum {dummy,RED, BLUE1, BLUE2,GREEN};


/* *** Constructors ********************************************************************************************/
Mesh::Mesh() : coordinates(), elements(), edge2nodes(), element2edges(), dirichlet(), neumann(),
               dirichlet2edges(), neumann2edges(), numEdges(0), numNodes(0), numElements(), numDirichlet(),
                     numNeumann(), coupling(), markedElements() 
{
    
}

Mesh::Mesh(Mesh &rhs)
{
        coordinates = rhs.coordinates;
        elements = rhs.elements; 
        dirichlet=rhs.dirichlet;
        neumann=rhs.neumann;
        edge2nodes = rhs.edge2nodes;
        element2edges = rhs.element2edges;
        dirichlet2edges = rhs.dirichlet2edges;
        neumann2edges=rhs.neumann2edges;
        coupling = rhs.coupling;
        markedElements = rhs.markedElements;
        numNodes = rhs.numNodes;
        numEdges = rhs.numEdges;
        numElements = rhs.numElements;
        numDirichlet = rhs.numDirichlet;
        numNeumann = rhs.numNeumann;
}


Mesh::Mesh(const Matrix &_coordinates, const IndexMatrix &_elements, const IndexMatrix &_dirichlet, 
               const IndexMatrix &_neumann, const IndexVector _elements2procs, const IndexMatrix _skeleton, 
               int _numCrossPoints_gl): 
           coordinates(), elements(), dirichlet(), neumann(), coupling(), markedElements()
{
    /* *** parallel version *** */
    if(_numCrossPoints_gl > 0){ 
      
        /* *** set local mesh for each processor *** */
        _distributeMesh(_coordinates, _elements , _dirichlet , _neumann,
                      _elements2procs, _skeleton, _numCrossPoints_gl);
                  
      /* *** set sizes *** */
        numNodes     = coordinates.numRows();    
        numElements  = elements.numRows();                                    
    numEdges     = numNodes+numElements-1; //Eulerscher Polyedersatz
        numDirichlet = dirichlet.size();
        numNeumann   = neumann.size();
    
        /* *** generate *2eges, edge2nodes (and local2globalCrossPoints) *** */
      provideGeometricDataMPI();

    /* *** serial version *** */
    }else{
      
      /* *** set coordinates elements dirichlet and neumann *** */
    coordinates = _coordinates;
    elements = _elements;
    _convertBdryData(_dirichlet, dirichlet);
    if(_neumann.numRows()>0){
      _convertBdryData(_neumann, neumann);
    }
    
      /* *** set sizes *** */
        numNodes     = coordinates.numRows();    
        numElements  = elements.numRows();                                    
    numEdges     = numNodes+numElements-1; //Eulerscher Polyedersatz
        numDirichlet = dirichlet.size();
        numNeumann   = neumann.size();
          
        /* *** generate *2eges, edge2nodes (and local2globalCrossPoints) *** */
      provideGeometricData();      
  }

}


Mesh::~Mesh(){
    
}


/* *** operators  *********************************************************************************************/

Mesh & Mesh::operator=(const Mesh &rhs)
{
  coordinates = rhs.coordinates;
    elements = rhs.elements; 
    dirichlet=rhs.dirichlet;
    neumann=rhs.neumann;
    edge2nodes = rhs.edge2nodes;
    element2edges = rhs.element2edges;
    dirichlet2edges = rhs.dirichlet2edges;
    neumann2edges=rhs.neumann2edges;
    coupling = rhs.coupling;
    markedElements = rhs.markedElements;
    numNodes = rhs.numNodes;
    numEdges = rhs.numEdges;
    numElements = rhs.numElements;
    numDirichlet = rhs.numDirichlet;
    numNeumann = rhs.numNeumann;
  
  return *this;
}


/* *** check if mesh is distributed (MPI) or not (serial) */
bool Mesh::distributed(){
    if (coupling.numCoupling==0) return false;
    else return true;
}

/* *** functions to create the data structure *2edges and edge2nodes (and local2globalCrossPoints) ************/

void Mesh::provideGeometricData()
{
    
    int nE = elements.numRows();
    int nD = 0;
    if(numDirichlet>0) for(int k=0;k<numDirichlet;k++) nD += dirichlet[k].length()-1;
    int nN = 0;
    if(numNeumann>0) for(int k=0;k<numNeumann;k++) nN += neumann[k].length()-1;
        
    /* *** build help Vectors 
    *         /  1st column elements  \                                        /  2nd column elements  \ 
    *         |  2nd column elements  |                                        |  3rd column elements  |
    *     I = |  3rd column elements  |       and     J = |  1st column elements  | 
    *         |  2nd column dirichlet |                                        |  1st column dirichlet |
    *         \  2nd column neumann  /                                        \  1st column neumann  /     
    *  every edges appears two times in I and J !
    */    
    IndexVector I(3*nE+nD+nN), J(3*nE+nD+nN);
    I.set(0      ,nE, elements.data()     );     /* 1st col */
    I.set(nE     ,nE, elements.data()+  nE);  /* 2nd col*/
    I.set(2*nE   ,nE, elements.data()+2*nE);     /* 3rd col */
    J.set(0      ,nE, elements.data()+  nE);     /* 2nd col */
    J.set(nE     ,nE, elements.data()+2*nE);  /* 3rd col */
    J.set(2*nE   ,nE, elements.data()     );  /* 1st col */
    
    // set Dirichlet
    if(numDirichlet>0){
        int index=0;
        for(int k=0; k<numDirichlet; k++){
            I.set(3*nE+index, dirichlet[k].length()-1, dirichlet[k].data() + 1);
            J.set(3*nE+index, dirichlet[k].length()-1, dirichlet[k].data()    );
            index += dirichlet[k].length()-1;
        }
    }    
    
    // set Neumann
    if(numNeumann>0){
        int index=0;
        for(int k=0; k<numNeumann; k++){
            I.set(3*nE+nD+index, neumann[k].length()-1, neumann[k].data() + 1);
            J.set(3*nE+nD+index, neumann[k].length()-1, neumann[k].data()    );
            index += neumann[k].length()-1;
        }
    }
    /* *** create IndexMatrix edge2nodes with 
           edge2nodes(i)=[j,k] <=> the idices of the nodes of the i-th edge are j and k*/
    int idxIJ=0,idxJI=0;
    IndexVector indexIJ1(numEdges),indexIJ2(numEdges);
    IndexVector indexJI1(numEdges),indexJI2(numEdges);
    IndexVector edgeNumber(I.length()); 
    IndexVector tmp1(numEdges), tmp2(numEdges);
    edge2nodes.resize(numEdges,2);
    for(int k=0;k<I.length();k++){
        if(I(k)<J(k)){
            /* *** extract edges with I<J */
            indexIJ1(idxIJ) = I(k)-1;
            indexIJ2(idxIJ) = J(k)-1;            
            tmp1(idxIJ) = idxIJ+1;
            /* *** fix numbering of edges with I<J  and write it in edgeNumber*/
            edgeNumber(k) = idxIJ+1;
            /* *** set edge2nodes */
            edge2nodes(idxIJ,0) = I(k);
            edge2nodes(idxIJ,1) = J(k);
            idxIJ++;
            
        }else{
            /* *** extract edges with J<I */
            indexJI1(idxJI) = J(k)-1;
            indexJI2(idxJI) = I(k)-1;
            tmp2(idxJI) = k+1;
            idxJI++;
        }
    }
    /* *** write number of edges with J<I in edgeNumber */ 
    CRSIndexMatrix numberingIJ(indexIJ1,indexIJ2 ,tmp1);
    CRSIndexMatrix idxJI2IJ(indexJI1,indexJI2,tmp2);
    for(int k=0;k<numEdges;k++){
        edgeNumber( idxJI2IJ(k)-1) = numberingIJ(k);
    }
    /* *** build element2edges*/
    element2edges.resize(nE,3);
    element2edges.set(0,nE,0,1,edgeNumber.data()      );
    element2edges.set(0,nE,1,1,edgeNumber.data()+   nE);
    element2edges.set(0,nE,2,1,edgeNumber.data()+2*nE );
    
    /* *** build dirichlet2edges and neumann2edges */
    dirichlet2edges.resize(numDirichlet),neumann2edges.resize(numNeumann);
    int index = 0;
    for(int k=0;k<numDirichlet;k++){
        int size=dirichlet[k].length()-1;
        dirichlet2edges[k].resize(size);
        dirichlet2edges[k].set(0, size, edgeNumber.data()+3*nE+index);
        index+= size;
    }
    index = 0;
    for(int k=0;k<numNeumann;k++){
        int size=neumann[k].length()-1;
        neumann2edges[k].resize(size);
        neumann2edges[k].set(0, size, edgeNumber.data()+3*nE+nD+index);
        index+= size;
    }
    
}

void Mesh::provideGeometricDataMPI()
{
    //--- Typedefs for FLENS
    // Used to set values from a C-Array
    typedef flens::DenseVector<flens::Array<int>::View> VecView;
    typedef flens::Array<int>::View View;
    //--- Typedefs for FLENS

    int nE = numElements;
    int nD = 0;
    if(numDirichlet>0) for(int k=0;k<numDirichlet;k++) nD += dirichlet[k].length()-1;
    int nN = 0;
    if(numNeumann>0) for(int k=0;k<numNeumann;k++) nN += neumann[k].length()-1;
    int nC = 0;  // number of edges on coupling boundaries
    for(int k=0;k<coupling.numCoupling;k++) nC+=coupling.boundaryNodes[k].length()-1;
    
    /* *** build help Vectors 
    *          / 1st column elements  \                                        /  2nd column elements  \ 
    *         |  2nd column elements  |                                        |  3rd column elements  |
    *     I = |  3rd column elements  |       and     J = |  1st column elements  | 
    *         |  2nd column dirichlet |                                        |  1st column dirichlet |
    *         |  2nd column neumann   |                                        |  1st column neumann   |
    *         \  2nd column coupling /                    \  1st column coupling /    
    *  every edges appears twice in I and J !
    */
    IndexVector I(3*nE+nD+nN+nC), J(3*nE+nD+nN+nC);
    I.set(0      ,nE, elements.data()     );     /* 1st col */
    I.set(nE     ,nE, elements.data()+  nE);  /* 2nd col*/
    I.set(2*nE   ,nE, elements.data()+2*nE);     /* 3rd col */
    J.set(0      ,nE, elements.data()+  nE);     /* 2nd col */
    J.set(nE     ,nE, elements.data()+2*nE);  /* 3rd col */
    J.set(2*nE   ,nE, elements.data()     );  /* 1st col */
    
    // set Dirichlet
    if(numDirichlet>0){
        int index=0;
        for(int k=0; k<numDirichlet; k++){
            I.set(3*nE+index, dirichlet[k].length()-1, dirichlet[k].data() + 1);
            J.set(3*nE+index, dirichlet[k].length()-1, dirichlet[k].data()    );
            index+=dirichlet[k].length()-1;
        }
    }
    // set Neumann
    if(numNeumann>0){
        int index=0;
        for(int k=0; k<numNeumann; k++){
            I.set(3*nE+nD+index, neumann[k].length()-1, neumann[k].data() + 1);
            J.set(3*nE+nD+index, neumann[k].length()-1, neumann[k].data()    );
            index+=neumann[k].length()-1;
        }
    }
    // set coupling
    if(coupling.numCoupling>0){
        int index=0;
        for(int k=0; k<coupling.numCoupling; k++){
            I.set(3*nE+nD+nN+index, coupling.boundaryNodes[k].length()-1, coupling.boundaryNodes[k].data() + 1);
            J.set(3*nE+nD+nN+index, coupling.boundaryNodes[k].length()-1, coupling.boundaryNodes[k].data()    );
            index+=coupling.boundaryNodes[k].length()-1;
        }
    }    

    /* *** create IndexMatrix edge2nodes with 
           edge2nodes(i)=[j,k] <=> the idices of the nodes of the i-th edge are j and k*/
    int idxIJ=0,idxJI=0;
    IndexVector indexIJ1(numEdges),indexIJ2(numEdges);
    IndexVector indexJI1(numEdges),indexJI2(numEdges);
    IndexVector edgeNumber(I.length()); 
    IndexVector tmp1(numEdges), tmp2(numEdges);
    edge2nodes.resize(numEdges,2);
    for(int k=0;k<I.length();k++){
        if( I(k)<J(k) ){
            /* *** extract edges with I<J */
            indexIJ1(idxIJ) = I(k)-1;
            indexIJ2(idxIJ) = J(k)-1;            
            tmp1(idxIJ) = idxIJ+1;
            /* *** fix numbering of edges with I<J  and write it in edgeNumber*/
            edgeNumber(k) = idxIJ+1;
            /* *** set edge2nodes */
            edge2nodes(idxIJ,0) = I(k);
            edge2nodes(idxIJ,1) = J(k);
            idxIJ++;
        }else{
            /* *** extract edges with J<I */
            indexJI1(idxJI) = J(k)-1;
            indexJI2(idxJI) = I(k)-1;
            tmp2(idxJI) = k+1;
            idxJI++;
        }
    }
    
    /* *** write number of edges with J<I in edgeNumber */     
    CRSIndexMatrix numberingIJ(indexIJ1,indexIJ2 ,tmp1);
    CRSIndexMatrix idxJI2IJ(indexJI1,indexJI2,tmp2);
    for(int k=0;k<numEdges;k++){
        edgeNumber( idxJI2IJ(k)-1) = numberingIJ(k);
    }
    /* *** build element2edges*/
    element2edges.resize(nE,3);
    element2edges.set(0,nE,0,1,edgeNumber.data()      );
    element2edges.set(0,nE,1,1,edgeNumber.data()+   nE);
    element2edges.set(0,nE,2,1,edgeNumber.data()+2*nE );
    /* *** build dirichlet2edges and neumann2edges */
    if(numDirichlet>0){
        dirichlet2edges.resize(numDirichlet);
        int index = 0;
        for(int k=0;k<numDirichlet;k++){
            int size=dirichlet[k].length()-1;
            dirichlet2edges[k].resize(size);
            dirichlet2edges[k].set(0, size, edgeNumber.data()+3*nE+index);
            index+= size;
        }
    }
    if(numNeumann>0){    
        int index = 0;
        neumann2edges.resize(numNeumann);
        for(int k=0;k<numNeumann;k++){
            int size=neumann[k].length()-1;
            neumann2edges[k].resize(size);
            neumann2edges[k].set(0, size, edgeNumber.data()+3*nE+nD+index);
            index+= size;
        }
    }    
    /* *** build coupling2edges */
    if(coupling.numCoupling>0){
        int index=0;
        coupling.coupling2edges.resize(coupling.boundaryNodes.size());
        for(int k=0;k<coupling.numCoupling;k++){
            int size = coupling.boundaryNodes[k].length()-1;
            coupling.coupling2edges[k].resize(size);
            coupling.coupling2edges[k] = VecView(View(size,edgeNumber.data()+3*nE+nD+nN+index));
            index+=size;
        }
    }
    
}

/* *** functions to write mesh ********************************************************************************/

void Mesh::writeData(int proc, std::string dir) 
{
  std::string strproc;
  if(proc==0){
      strproc="";
  }
  else{
      std::stringstream ss;
      ss << proc;
      strproc=ss.str();
  }

    std::string Coordinates= dir + "coordinates" + strproc + ".dat";
    coordinates.write(Coordinates);

    std::string Elements=dir + "elements" + strproc + ".dat";
    elements.write(Elements);

    std::fstream f;
    std::string Dirichlet = dir + "dirichlet"+ strproc + ".dat";
    f.open(Dirichlet.c_str(), std::ios::out);
    for(int k=0;k<numDirichlet;k++){
        for(int j=0;j<dirichlet[k].length()-1;j++)
            f << dirichlet[k](j)<<"   "<<dirichlet[k](j+1) << std::endl;
    }
    f.close();

    std::string Neumann=dir + "neumann"+ strproc + ".dat";
    f.open(Neumann.c_str(), std::ios::out);
    for(int k=0;k<numNeumann;k++){
        for(int j=0;j<neumann[k].length()-1;j++)
            f << neumann[k](j)<<"   "<< neumann[k](j+1) << std::endl;
    }
    f.close();    

    std::string Edge2nodes=dir + "edge2nodes" + strproc + ".dat";
    edge2nodes.write(Edge2nodes);

    std::string Element2edges=dir + "element2edges" + strproc + ".dat";
    element2edges.write(Element2edges);

    std::string Dirichlet2edges=dir + "dirichlet2edges" + strproc + ".dat";
    f.open(Dirichlet2edges.c_str(), std::ios::out);
    for(int k=0;k<numDirichlet;k++){
        for(int j=0;j<dirichlet2edges[k].length();j++)
            f << dirichlet2edges[k](j) << std::endl;
    }
    f.close();

    std::string Neumann2edges=dir + "neumann2edges" + strproc + ".dat";
    f.open(Neumann2edges.c_str(), std::ios::out);
    for(int k=0;k<numNeumann;k++){
        for(int j=0;j<neumann2edges[k].length();j++)
            f << neumann2edges[k](j) << std::endl;
    }
    f.close();

}

/* *** functions to refine the mesh ***************************************************************************/

/* red refinement: each element is split into four elements */
void Mesh::refineRed()
{
        
  /* *** compute new coordinates  */
  Matrix oldCoordinates(coordinates);
  coordinates.resize(numNodes+numEdges,2);
  coordinates.set(0,numNodes,0, 2, oldCoordinates.data()) ;
  for (int j=0; j<numEdges; j++) { // new Coordinates
    // set x-coordinates
    coordinates(numNodes+j,0) = (1./2.)*( coordinates( edge2nodes(j,0)-1, 0) +
                                          coordinates( edge2nodes(j,1)-1, 0)  );
    // set y-coordinates
    coordinates(numNodes+j,1) = (1./2.)*( coordinates( edge2nodes(j,0)-1, 1) +
                                          coordinates( edge2nodes(j,1)-1, 1)  );
  }
    /* *** compute new elements  */
  IndexMatrix oldElements(elements);
  elements.resize(numElements*4,3);

  for (int j=0; j<numElements; j++) { // elementwise refinement
    // first new element
    elements(4*j,0) = oldElements(j,0);
    elements(4*j,1) = numNodes+element2edges(j,0);
    elements(4*j,2) = numNodes+element2edges(j,2);

    //  second new element
    elements(1 + 4*j,0) =   numNodes+element2edges(j,0);
    elements(1 + 4*j,1) =   oldElements(j,1);
    elements(1 + 4*j,2) =   numNodes+element2edges(j,1);

    //  third new element
    elements(2 + 4*j,0) =  numNodes+element2edges(j,1);
    elements(2 + 4*j,1) =  oldElements(j,2);
    elements(2 + 4*j,2) =  numNodes+element2edges(j,2);

    //  fourth new element
    elements(3 + 4*j,0) =   numNodes+element2edges(j,2);
    elements(3 + 4*j,1) =   numNodes+element2edges(j,0);
    elements(3 + 4*j,2) =   numNodes+element2edges(j,1);
  }
    
    /* *** set new dirichlet boundary */
    if(numDirichlet>0){
        std::vector<IndexVector> oldDirichlet(dirichlet);
        for(int k=0;k<numDirichlet;k++){
            dirichlet[k].resize(2*oldDirichlet[k].length()-1);
              for (int j=0; j<oldDirichlet[k].length()-1; j++) {
            dirichlet[k](2*j)   =  oldDirichlet[k](j);
            dirichlet[k](2*j+1) =  dirichlet2edges[k](j)+numNodes;
              }
            dirichlet[k](dirichlet[k].length()-1) = oldDirichlet[k]( oldDirichlet[k].length()-1 );
        }
    }

    /* *** set new neumann boundary */
    if(numNeumann>0){
        std::vector<IndexVector> oldNeumann(neumann);
        for(int k=0;k<numNeumann;k++){
            neumann[k].resize(2*oldNeumann[k].length()-1);
          for (int j=0; j<oldNeumann[k].length()-1; j++) {
          neumann[k](2*j)   =  oldNeumann[k](j);
          neumann[k](2*j+1) =  neumann2edges[k](j)+numNodes;
          }
            neumann[k]( neumann[k].length()-1 ) = oldNeumann[k]( oldNeumann[k].length()-1 );
        }
    }
    
    /* *** set new coupling boundary => only the field boundaryNodes has to be adjusted */
    if (coupling.numCoupling>0) {
       Coupling oldCoupling(coupling);

      for (int k=0; k<coupling.numCoupling; k++){    //for each coupling boundary
          int j; 
          coupling.boundaryNodes[k].resize(2*oldCoupling.boundaryNodes[k].length()-1); 
          for (j=0; j<oldCoupling.boundaryNodes[k].length()-1; j++){
        coupling.boundaryNodes[k](2*j+1)    =  oldCoupling.boundaryNodes[k](j+1);
          coupling.boundaryNodes[k](2*j+2)  =  coupling.coupling2edges[k](j+1)+numNodes;
      }
          coupling.boundaryNodes[k]( coupling.boundaryNodes[k].length() ) = 
            oldCoupling.boundaryNodes[k]( oldCoupling.boundaryNodes[k].length() );                                                            
    }
  }

    /* *** set new sizes */
    numNodes     = coordinates.numRows();
    numElements  = elements.numRows();
    numEdges     = numNodes+numElements-1; //Eulerscher Polyedersatz
    
    // provide geometric data => *2edges and edge2nodes
    if(coupling.numCoupling==0) provideGeometricData();
    else provideGeometricDataMPI();

}

/* adaptive rgb refinement */
IndexVector Mesh::refineRGB(Vector &material)
{ 
    IndexVector dummy(1);
    return dummy;
}


/* *** private functions *************************************************************************************/

/* distribute the mesh such that each processor has its local mesh */
/* each process extracts its local mesh from the global one */
void Mesh::_distributeMesh(const Matrix &coordinates_gl, const IndexMatrix &elements_gl, const IndexMatrix &dirichlet_gl,
                                     const IndexMatrix &neumann_gl, const IndexVector &elements2procs, const IndexMatrix &skeleton,
                                     int numCrossPoints_gl)
{
    
    const int rank     = MPI::COMM_WORLD.Get_rank();
  const int numProcs = MPI::COMM_WORLD.Get_size();

    /* *** Compute local elements for each processor  *************************************************/
    // find indices of the elements of processor rank; 
    int numElem=0;
    for(int k=0;k<elements2procs.length();k++)
        if(elements2procs(k)==(rank+1) ) numElem++;
        
    IndexVector localInd(numElem);    
    int pos=0;
    for(int k=0;k<elements2procs.length();k++)
        if(elements2procs(k)==(rank+1)){
            localInd(pos) = k;
            pos++;
    }
  // Elements_loc covers the local elements but in global(!) index.
  int numElements_loc = localInd.length();
  elements.resize(numElements_loc, 3);
  for (int i=0; i<numElements_loc; ++i) {
        elements(i,0) = elements_gl( localInd(i) , 0);
        elements(i,1) = elements_gl( localInd(i) , 1);
        elements(i,2) = elements_gl( localInd(i) , 2);
  }
  
  /* *** Compute mapping from global to local indexing and set local2globalCrossPoints *************/
    IndexVector global2local(coordinates_gl.numRows());
  // First, we set cross points and compute number of neighbouring processes
  int locIndex = 0;
    int numNeighbours = 0;
    IndexVector colors(numProcs);
  for (int i=0; i<skeleton.numRows(); i++) {
      if ( skeleton(i,2)==(rank+1) || skeleton(i,3) == (rank+1) ) {
                    // here, we use the fact that the cross points are listed first!
          if (global2local( skeleton(i,0)-1 )==0  && skeleton(i,0) <= numCrossPoints_gl) {
              global2local( skeleton(i,0)-1 ) = ++locIndex;    
          }
          if (global2local( skeleton(i,1)-1 )==0  && skeleton(i,1) <= numCrossPoints_gl) {
              global2local( skeleton(i,1)-1 ) = ++locIndex;
          }
                    // build Vector with colors, so that every color only occurs once
                    bool newcolor=true;
                    for(int k=0;k<numNeighbours;k++)
                        if ( colors(k)==skeleton(i,4) ) newcolor=false;
                    
                    if(newcolor){
                        colors(numNeighbours) = skeleton(i,4);
                        numNeighbours++;
                    }
      }
  }
    int numCrossPoints_loc = locIndex;
    
    // set other nodes in global2local
    for (int i=0; i<numElements_loc; i++) {
      for (int j=0; j<3; j++){
          if (global2local( elements(i,j)-1 ) == 0) {
              global2local( elements(i,j)-1 ) = ++locIndex;
          }
      }
  }

  int numCoordinates_loc = locIndex;
    /* *** generate local skeleton *********************************************************************/
    IndexMatrix skeleton_loc(skeleton.numRows(),5);
    for (int i=0; i<skeleton.numRows(); i++) {
        if ( global2local( skeleton(i,0)-1 )!=0 && global2local( skeleton(i,1)-1 )!=0 ){ // local edge
            skeleton_loc(i,0) = global2local( skeleton(i,0)-1 );
            skeleton_loc(i,1) = global2local( skeleton(i,1)-1 );
            skeleton_loc(i,2) = skeleton(i,2);
            skeleton_loc(i,3) = skeleton(i,3);
            skeleton_loc(i,4) = skeleton(i,4);
        }
    }
    
    /* *** set coupling  *******************************************************************************/ 
    // initialize coupling
    coupling.numCrossPoints = numCrossPoints_gl;
    coupling.numCoupling = numNeighbours;
    coupling.boundaryNodes.resize(numNeighbours);
    coupling.local2globalCrossPoints.resize(numCrossPoints_loc);
    coupling.crossPointsBdryData.resize(numCrossPoints_loc);
    coupling.crossPointsNumProcs.resize(numCrossPoints_loc);
    coupling.colors.resize(numNeighbours);
    coupling.neighbourProcs.resize(numNeighbours);
    /* *** set crossPointsBdryData  ********************************************************************/ 
    // here, we use the fact that cross points are listed first in local indexing
    for(int j=0;j<dirichlet_gl.numRows(); j++){
            if( dirichlet_gl(j,0) <= numCrossPoints_gl && global2local( dirichlet_gl(j,0)-1 ) != 0 ){ 
                coupling.crossPointsBdryData( global2local( dirichlet_gl(j,0)-1 )) = 1;
            }    
                
            if( dirichlet_gl(j,1) <= numCrossPoints_gl && global2local( dirichlet_gl(j,1)-1 ) != 0 ){
                coupling.crossPointsBdryData( global2local( dirichlet_gl(j,1)-1 )) = 1;    
            }    
    }
    
    /* *** set boundary nodes ( boundary nodes have to be ordered!) and local2globalCrossPoints ********/
    // we start at one cross point and walk through skeleton_loc to next cross point 
    // => indices for boundary nodes 
    IndexVector node2Skeleton(numCoordinates_loc);
    IndexVector ptr(skeleton.numRows()+1);
    
    // initialize node2Skeleton
    for(int k=0;k<skeleton_loc.numRows();k++){
        if(skeleton_loc(k,0)!=0){
            node2Skeleton( skeleton_loc(k,0)-1 ) = k+1;
        }    
    }
    int index=0, sum=1, numBoundaryNodes=0;
    while(sum!=0) { 
        // initialize  ptr, we start with the row of first cross point, that we find in skeleton_loc
        for(int k=0;k<skeleton_loc.numRows();k++){
            if (skeleton_loc(k,0)!=0 && skeleton_loc(k,0) <= numCrossPoints_loc){
                ptr(0) = k+1;
                //set local2globalCrossPoints
                coupling.local2globalCrossPoints( skeleton_loc(k,0) ) = skeleton(k,0);
                break;
            }
        }
        // walk through skeleton until the next cross points and save indices of the rows in ptr                     
        for(int k=1;k<=skeleton.numRows();k++){
            if (skeleton_loc( ptr(k-1)-1,1 ) <= numCrossPoints_loc){ // we reached the next crosspoint! 
                numBoundaryNodes=k+1;
                //set local2globalCrossPoints, if it is not already set
                if( coupling.local2globalCrossPoints( skeleton_loc( ptr(k-1)-1,1 ) ) == 0 )
                      coupling.local2globalCrossPoints( skeleton_loc( ptr(k-1)-1,1 ) ) =  skeleton( ptr(k-1)-1,1 );
                
                break;
            }
            ptr(k) = node2Skeleton( skeleton_loc( ptr(k-1)-1,1 )-1);
        }
        // set color    
        coupling.colors(index+1) = skeleton_loc( ptr(0)-1,4 );
        if ( skeleton_loc(ptr(0)-1,2)==(rank+1) ){    // skeleton elements are ordered counter clockwise    
            // set neighbourProcs
            coupling.neighbourProcs(index+1) = skeleton_loc( ptr(0)-1,3);
            //set boundary nodes ( and delete skeleton)
            coupling.boundaryNodes[index].resize(numBoundaryNodes);
            
            for(int k=0;k<numBoundaryNodes-1;k++){
                coupling.boundaryNodes[index](k+1) = skeleton_loc( ptr(k)-1,0);
                if(k==numBoundaryNodes-2) 
                    coupling.boundaryNodes[index](k+2) = skeleton_loc( ptr(k)-1,1);
                //set skeleton to zero
                skeleton_loc( ptr(k)-1, 0) = 0;
                skeleton_loc( ptr(k)-1, 1) = 0;
            }            
        }    else if ( skeleton_loc(ptr(0)-1,3)==(rank+1) ){  // skeleton elements are ordered clockwise => switch    
            // set neighbourProcs
            coupling.neighbourProcs(index+1) = skeleton_loc( ptr(0)-1,2 );
            // set boundary nodes ( and delete skeleton)
            coupling.boundaryNodes[index].resize(numBoundaryNodes);
            for(int k=0;k<numBoundaryNodes-1;k++){ 
                coupling.boundaryNodes[index](k+1) = skeleton_loc( ptr(numBoundaryNodes-k-2)-1,1 );
                if(k==numBoundaryNodes-2)
                    coupling.boundaryNodes[index](k+2) = skeleton_loc( ptr(0)-1,0 );
                //set skeleton to zero
                skeleton_loc( ptr(numBoundaryNodes-k-2)-1, 0) = 0;
                skeleton_loc( ptr(numBoundaryNodes-k-2)-1, 1) = 0;
            }
        }
        // compute sum of first column (if 0 all coupling boundarys are done!)
        sum=0;
        for(int k=0;k<skeleton_loc.numRows();k++)
            sum += skeleton_loc(k,0);
        
        index++;    
    }
    
    /* *** build local coordinates, elements *******************************************************************/
    // elements
    for (int i=0; i<numElements_loc; i++) {
      elements(i,0) = global2local( elements(i,0)-1 );
      elements(i,1) = global2local( elements(i,1)-1 );
        elements(i,2) = global2local( elements(i,2)-1 );
    }
    // coordinates
    coordinates.resize(numCoordinates_loc,2);
    for (int i=0; i<global2local.length(); i++) {
    if (global2local(i)!=0) {
        coordinates( global2local(i)-1,0) = coordinates_gl(i,0);
        coordinates( global2local(i)-1,1) = coordinates_gl(i,1);
    }
    }


    /* *** build local dirichlet/neumann (a vector that contains the SORTED Dirichlet/Neumann nodes ) ******/
    //convert Dirichlet from global to local indexing
    IndexMatrix tmpdir(dirichlet_gl); 
    for(int k=0;k<dirichlet_gl.numRows();k++){
        if( global2local( tmpdir(k,0)-1 )!=0 && global2local( tmpdir(k,1)-1) !=0){
            dirichlet_gl(k,0) = global2local( tmpdir(k,0)-1 );
            dirichlet_gl(k,1) = global2local( tmpdir(k,1)-1 );
        }else{
            dirichlet_gl(k,0) = 0;
            dirichlet_gl(k,1) = 0;    
        }    
    }
  _convertBdryData(dirichlet_gl, dirichlet);
  
    // no Neumann cycle possible since we do not calculate Neumann problems, i.e. |Gamma_D|>0 ! 
    if(neumann_gl.numRows()>0){
        
        //convert Neumann from global to local indexing
        IndexMatrix tmpN(neumann_gl);
        for(int k=0;k<neumann_gl.numRows();k++){
            if( global2local( tmpN(k,0)-1 ) != 0 && global2local( tmpN(k,1)-1 ) != 0 ){
                neumann_gl(k,0) = global2local(tmpN(k,0)-1);
                neumann_gl(k,1) = global2local(tmpN(k,1)-1);
            }else{
                neumann_gl(k,0) = 0;
                neumann_gl(k,1) = 0;    
            }    
        }
        _convertBdryData(neumann_gl, neumann);
    }    
        
     
    
    
    /* *** build maxColor ***********************************************/
    int maxColorLoc = 0;
    for (int i=1; i<=coupling.colors.length(); ++i)
    {
        if(coupling.colors(i)>maxColorLoc) maxColorLoc=coupling.colors(i);
    }

    MPI::COMM_WORLD.Allreduce(&maxColorLoc, &coupling.maxColor, 1, MPI::INT, MPI::MAX);
    /* *** build coupling_loc.crossPointsNumProcs **********************************************************/
    for (int j=1 ; j<=numCrossPoints_loc; j++) {
          for (int k=0; k<skeleton.numRows(); k++) {
              if (skeleton(k,0) == coupling.local2globalCrossPoints(j) ||
                  skeleton(k,1) == coupling.local2globalCrossPoints(j) )
                            coupling.crossPointsNumProcs(j) +=1;
          }
      }
    for (int j=1 ; j<=numCrossPoints_loc; j++)
        if(coupling.crossPointsNumProcs(j)==1)
            coupling.crossPointsNumProcs(j) = 2;

}

/* convert dirichlet and neumann data from n x 2 matrices to vectors of IndexVectors */
/* each connected dirichlet or neumann boundary is one index vector */
void Mesh::_convertBdryData(const IndexMatrix bdry_tmp, std::vector<IndexVector> &bdry)
{

  int numCoordinates = coordinates.numRows();
  IndexVector ptrD1(bdry_tmp.numRows()+1), ptrD2(bdry_tmp.numRows()+1);
  // initialize node2Bdry1 and node2Bdry2
  IndexVector node2Bdry1(numCoordinates );
  IndexVector node2Bdry2(numCoordinates );
  for(int k=0;k<bdry_tmp.numRows();k++){
      if(bdry_tmp(k,0)!=0){
          node2Bdry1(bdry_tmp(k,0)-1) = k+1;
          node2Bdry2(bdry_tmp(k,1)-1) = k+1;
      }    
  }
  //compute sum (check if there are local bdry edges)
  int sumD=0;
  for(int k=0;k<bdry_tmp.numRows();k++){
      sumD+= bdry_tmp(k,0);
  }
  while(sumD!=0) { 
      // initialize  ptrD1, we start with the row of first local Bdry node
      for(int k=0;k< bdry_tmp.numRows();k++){
          if (bdry_tmp(k,0)!=0){
              ptrD1(0) = k+1;
              break;
          }
      }
      int numD1=0, numD2=0;
      bool foundCycle=false;
      // walk through Bdry until we reach start node                      
      for(int k=1;k<=bdry_tmp.numRows();k++){
          ptrD1(k) = node2Bdry1(bdry_tmp(ptrD1(k-1)-1 ,1) -1 );
          numD1++;
          // we reached start point again => finished
          if (ptrD1(k) == ptrD1(0) ){  
              foundCycle=true;
              break;
          }
          // we don't have a connected piece of boundary => walk other direction
          if (ptrD1(k) == 0 ) break;
      }
      // walk through Bdry in other direction
      if(foundCycle==false){
          // initialize  ptrD2
          ptrD2(0) = ptrD1(0);
          for(int k=1;k<=bdry_tmp.numRows();k++){
              numD2++;
              ptrD2(k) = node2Bdry2(    bdry_tmp(ptrD2(k-1)-1 ,0)-1 );
              if (ptrD2(k) == 0){  
              break;
              }
          }
          numD2--; // we have first entry twice    
      }
      // save bdry nodes in bdrytmp and set bdry_tmp to zero
      IndexVector bdry2(numD1+numD2+1);
      for(int k=0; k<numD2; k++){
          bdry2(k) =  bdry_tmp(ptrD2(numD2-k)-1,0);
          bdry_tmp(ptrD2(numD2-k)-1, 0) = 0;
      }
      for(int k=0; k<numD1; k++){
          bdry2(k+numD2)  = bdry_tmp(ptrD1(k)-1,0);
          bdry_tmp( ptrD1(k)-1, 0) = 0;
      }
    
      bdry2(numD1+numD2) = bdry_tmp(ptrD1(numD1-1)-1,1);
    
      // resize bdry_loc
      int size = bdry.size();
      std::vector<IndexVector> bdry_old(bdry);
      
      bdry.resize(size+1);
      
      for(int k=0;k<size;k++) bdry[k]= bdry_old[k];
      // set Bdry
      bdry[size] = bdry2;
      // compute sum
      sumD=0;
      for(int k=0;k<bdry_tmp.numRows();k++){
          sumD+= bdry_tmp(k,0);
      }    
  }

}

/* *** Functions to read the mesh from .dat-files  (serial and parallel) *************************************/
void readMesh(char * directory,    Matrix &coordinates, IndexMatrix &elements, IndexMatrix & dirichlet,
              IndexMatrix &neumann)
{
            
    char tmp[256];
    strcpy(tmp,directory);
    coordinates.read(strcat(tmp,"/coordinates.dat"));
    strcpy(tmp,directory);
    elements.read(strcat(tmp,"/elements.dat"));
    strcpy(tmp,directory);
    dirichlet.read(strcat(tmp,"/dirichlet.dat"));
    strcpy(tmp,directory);
    neumann.read(strcat(tmp,"/neumann.dat"));        
}

void readMeshMPI(char * input, Matrix &coordinates,IndexMatrix &elements, IndexMatrix & dirichlet,
                     IndexMatrix &neumann,IndexVector &elements2procs, IndexMatrix &skeleton, int &numCrossPoints_gl)
{
        
    IndexMatrix elementsold,dirichletold,neumannold, skeletonold;
    Matrix coordinatesold;

    char directory[256] = "";
    char input_file[256] ="";
    strcat(directory,input);
  /* **** processor 0 reads the geometry and distributes it     */
  /* **** to all other processors.                              */
  /* initialize mesh */
  strcpy(input_file,directory);
    coordinatesold.read(strcat(input_file,"/coordinates.dat"));
    strcpy(input_file,directory);
    elementsold.read(strcat(input_file,"/elements.dat"));
    strcpy(input_file,directory);
    dirichletold.read(strcat(input_file,"/dirichlet.dat"));
    strcpy(input_file,directory);
    neumannold.read(strcat(input_file,"/neumann.dat"));
    strcpy(input_file,directory);
    skeletonold.read(strcat(input_file,"/skeleton.dat"));
    strcpy(input_file,directory);
    elements2procs.read(strcat(input_file,"/elements2procs.dat"));

    /* ***rearrange mesh s.t. cross points are listed first ************************************************** */
    // build vector old2new to convert from old order to new order 
    int *old2new   = new int[coordinatesold.numRows()];
    memset(old2new, 0, coordinatesold.numRows()*sizeof(int));
    for(int k=0; k<skeletonold.numRows();k++){
        old2new[skeletonold(k,0)-1]++;
        old2new[skeletonold(k,1)-1]++;
    }
    for(int k=0;k<dirichletold.numRows();k++){
        old2new[dirichletold(k,0)-1]++;
        old2new[dirichletold(k,1)-1]++;
    }
    for(int k=0;k<neumannold.numRows();k++){
        old2new[neumannold(k,0)-1]++;
        old2new[neumannold(k,1)-1]++;
    }
    numCrossPoints_gl=0;
    for(int k=0;k<coordinatesold.numRows();k++)
        if(old2new[k]!=2 && old2new[k]!=0){ // cross point found!
            numCrossPoints_gl++;
        }
    int indexCrossPoints=1;
    int indexOtherPoints=numCrossPoints_gl+1;        
    for(int k=0;k<coordinatesold.numRows();k++){
        if(old2new[k]!=2 && old2new[k]!=0)     old2new[k]=indexCrossPoints++;
        else                                 old2new[k]=indexOtherPoints++;
    }    
    // set new dirichlet
    dirichlet.resize(dirichletold.numRows(),2);
    for(int k=0;k<dirichletold.numRows();k++){
        dirichlet(k,0) = old2new[dirichletold(k,0)-1];
        dirichlet(k,1) = old2new[dirichletold(k,1)-1];
    }
    // set new neumann
    neumann.resize(neumannold.numRows(),2);
    for(int k=0;k<neumannold.numRows();k++){
        neumann(k,0) = old2new[neumannold(k,0)-1];
        neumann(k,1) = old2new[neumannold(k,1)-1];
    }
    // set new elements
    elements.resize(elementsold.numRows(),3);
    for(int k=0;k<elementsold.numRows();k++){
        elements(k,0) = old2new[elementsold(k,0)-1];
        elements(k,1) = old2new[elementsold(k,1)-1];
        elements(k,2) = old2new[elementsold(k,2)-1];
    }
    // set new skeleton
    skeleton.resize(skeletonold.numRows(),5);
    for(int k=0;k<skeletonold.numRows();k++){
        skeleton(k,0) = old2new[skeletonold(k,0)-1];
        skeleton(k,1) = old2new[skeletonold(k,1)-1];
    }
    // set other columns of skeleton
    skeleton.set(0,skeletonold.numRows(),2,3,skeletonold.data()+ 2*skeletonold.numRows());
    // set new coordinates
         coordinates.resize(coordinatesold.numRows(),2);
    for(int k=0;k<coordinatesold.numRows();k++){
        coordinates(old2new[k]-1,0) = coordinatesold(k,0);
        coordinates(old2new[k]-1,1) = coordinatesold(k,1);
    }

}





