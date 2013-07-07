#ifndef GS_NOMPI_BLAS_CPP
#define GS_NOMPI_BLAS_CPP 1

#include "gs_nompi_blas.h"

//FLENS-based GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_dense_nompi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
   int    maxIterations = std::numeric_limits<int>::max(),
   double tol = std::numeric_limits<double>::epsilon())
{
  using namespace flens;

	typedef typename VB::ElementType  ElementType;
  typedef typename VB::IndexType    IndexType;
  typedef typename VB::NoView       VectorType;

  ElementType rNormSquare, tmp;
  VectorType r, bc_node(x.length());
  const ElementType  Zero(0), One(1);

  // Compute error norm
  blas::copy(b, r);
  blas::mv(NoTrans, -One, A, x, One, r);
  // Set residual to zero at dirichlet nodes
  for (int i=0; i<bc.length(); ++i)
  {
      r(bc(i)) = Zero;
      bc_node(bc(i)) = One;
  }
  rNormSquare = blas::dot(r, r);

  for (int k=1; k<=maxIterations; ++k)
  {
      // STOP criteria
      if (sqrt(rNormSquare)<=tol) return k-1;

      // Update vector
      for (IndexType j=x.firstIndex(), i=A.firstRow(); j<=x.lastIndex(); ++j, ++i)
      {
          // Update only free nodes
          if (bc_node(j)==Zero)
          {
              tmp = Zero;
              for (IndexType lA=A.firstCol(), lx=x.firstIndex(); lA<=A.lastCol(); ++lA, ++lx)
              {
                  tmp += A(i,lA)*x(lx);
              }
              x(j) += (b(j)-tmp)/A(i,i);
          } 
      }

      // Compute error norm
      blas::copy(b, r);
      blas::mv(NoTrans, -One, A, x, One, r);

      // Set residual to zero at dirichlet nodes
      for (int i=0; i<bc.length(); ++i)
      {
          r(bc(i)) = Zero;
      }
      rNormSquare = blas::dot(r, r);
  }

  // Max iterations reached
	return maxIterations;
};


//Wrapper: Funken --> FLENS --> Funken
int
gs_dense_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
							int maxIt, double tol)
{
  typedef int                                          IndexType;
  typedef flens::IndexOptions<IndexType, 1>            IndexBase;
  typedef flens::DenseVector<flens::Array<double> >    DenseVector;
    
  //Check if sizes of matrices & vectors fit:
  assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
 



 //##### This part has to be changed for dense matrices!!!!!!!!!!!!!!!!   
  //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
  flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
  funk2flens_CRSmat(fk_A, fl_A);

  // Densify
  flens::GeMatrix<flens::FullStorage<double> > Fl_A = fl_A;
 //##### This part has to be changed for dense matrices!!!!!!!!!!!!!!!!

 

 

  //Convert Funken Vector b --> FLENS DenseVector:
  DenseVector fl_b(fk_b.length());
  funk2flens_Vector(fk_b, fl_b);
    
  //Solve using the FLENS-based GS solver:
  int iterCount;
  DenseVector fl_x(fl_b.length());
  iterCount = gs_dense_nompi_blas(Fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

  //Convert solution FLENS DenseVector x --> Funken Vector:
  flens2funk_Vector(fl_x, fk_x);
  
  return iterCount;
};

#endif	//GS_NOMPI_BLAS_CPP