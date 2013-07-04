#ifndef SOLVER_H_
#define SOLVER_H_

#include <time.h>
#include <sys/time.h>
#include <vector>

#include "../LinearAlgebra/LinAlgHeader.hpp"

#include "MathOperationsMPI.hpp"
#include "Mesh.hpp"
#include "DataVector.hpp"
#include "Coupling.hpp"

/* *** CG-methods **************************************************************************************/
int CG(CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes,  int maxIt, double tol);

int CG_MPI(CRSMatrix &A, DataVector &x, DataVector &b, IndexVector &dirichletNodes,  int maxIt, double tol);
	
/* *** Gauß-Seidel-methods **************************************************************************************/
int forwardGS( CRSMatrix &A, Vector &x, Vector &b, IndexVector &dirichletNodes, int maxIt, double tol);

int forwardGS_MPI( CRSMatrix &A, DataVector &x, DataVector &b,Coupling &coupling,
                   IndexVector &dirichletNodes, int maxIt);
			
void solveTridiag(Vector &ldiag, Vector &diag, Vector &udiag, Vector &x, Vector &b);

#endif // SOLVER_H

