#ifndef GS_NOMPI_BLAS_CPP
#define GS_NOMPI_BLAS_CPP 1

#include "gs_nompi_blas.h"

<<<<<<< HEAD
//FLENS-based dense GS solver:
=======

//FLENS-based GS solver:
>>>>>>> chris
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
      for (IndexType i=1; i<=x.lastIndex(); ++i)
      {
          // Update only free nodes
          if (bc_node(i)==Zero)
          {
              tmp = Zero;
              for (IndexType j=1; j<=A.lastCol(); ++j)
              {
                  tmp += A(i,j)*x(j);
              }
              x(i) += (b(i)-tmp)/A(i,i);
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


//FLENS-based sparse GS solver:
template <typename MA, typename VX, typename VB, typename VBC>
int
gs_nompi_blas(const MA &A, const VB &b, VX &x, VBC &bc,
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

  // Access values of A
  const auto &_rows = A.engine().rows();
  const auto &_cols = A.engine().cols();
  const auto &_vals = A.engine().values();

  // Iteration start
  for (int k=1; k<=maxIterations; ++k)
  {
      // STOP criteria
      if (sqrt(rNormSquare)<=tol) return k-1;

      // Update vector
      for (IndexType i=1; i<=x.lastIndex(); ++i)
      {
          // Update only free nodes
          if (bc_node(i)==Zero)
          {
              tmp = Zero;
              ElementType Aii = Zero;
              for (IndexType j=_rows(i); j<_rows(i+1); ++j)
              {
                  tmp += _vals(j)*x(_cols(j));
                  if (_cols(j)==i) {
                      Aii = _vals(j);
                  }
              }
              assert(Aii!=Zero);
              x(i) += (b(i)-tmp)/Aii;
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
template <typename MA, typename VX, typename VB, typename VBC>
int
<<<<<<< HEAD
gs_dense_nompi_blas_wrapper(Matrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
=======
gs_dense_nompi_blas_wrapper(MA &fk_A, VX &fk_x, VB &fk_b, VBC &fk_bc,
>>>>>>> chris
							int maxIt, double tol)
{
  typedef int                                          IndexType;
  typedef flens::IndexOptions<IndexType, 1>            IndexBase;
  typedef flens::DenseVector<flens::Array<double> >    DenseVector;
    
  //Check if sizes of matrices & vectors fit:
  assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
   
  //Convert Funken Matrix A --> FLENS GeMatrix:
  flens::GeMatrix<flens::FullStorage<double> > fl_A(fk_A.numRows(),fk_A.numCols());
  funk2flens_mat(fk_A, fl_A);

  //Convert Funken Vector b --> FLENS DenseVector:
  DenseVector fl_b(fk_b.length());
  funk2flens_Vector(fk_b, fl_b);
    
  //Solve using the FLENS-based GS solver:
  int iterCount;
  DenseVector fl_x(fl_b.length());
  iterCount = gs_dense_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

  //Convert solution FLENS DenseVector x --> Funken Vector:
  flens2funk_Vector(fl_x, fk_x);
  
  return iterCount;
};


//Wrapper: Funken --> FLENS --> Funken
int
gs_nompi_blas_wrapper(CRSMatrix &fk_A, Vector &fk_x, Vector &fk_b, IndexVector &fk_bc,
              int maxIt, double tol)
{
  typedef int                                          IndexType;
  typedef flens::IndexOptions<IndexType, 1>            IndexBase;
  typedef flens::DenseVector<flens::Array<double> >    DenseVector;
    
  //Check if sizes of matrices & vectors fit:
  assert(fk_A.numRows()==fk_b.length() && fk_A.numCols()==fk_x.length());
 
  //Convert Funken CRSMatrix A --> FLENS CRS Matrix:
  flens::GeCRSMatrix<flens::CRS<double, IndexBase> > fl_A;
  funk2flens_CRSmat(fk_A, fl_A);

  //Convert Funken Vector b --> FLENS DenseVector:
  DenseVector fl_b(fk_b.length());
  funk2flens_Vector(fk_b, fl_b);
    
  //Solve using the FLENS-based GS solver:
  int iterCount;
  DenseVector fl_x(fl_b.length());
  iterCount = gs_nompi_blas(fl_A, fl_b, fl_x, fk_bc, maxIt, tol);

  //Convert solution FLENS DenseVector x --> Funken Vector:
  flens2funk_Vector(fl_x, fk_x);
  
  return iterCount;
};

#endif	//GS_NOMPI_BLAS_CPP