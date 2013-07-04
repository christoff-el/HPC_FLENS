#ifndef MATHOPERATIONSMPI_H_
#define MATHOPERATIONSMPI_H_

#include <math.h>
#include <mpi.h>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "DataVector.hpp"
#include "Mesh.hpp"


/* *** Vector operations ********************************************************/
double max(DataVector& u);
double norm(DataVector& u);

/* *** Vector-Vector operations *************************************************/
double dot(DataVector &uTypeI, DataVector &vTypeII);

// u += alpha*v
void add(DataVector &u, DataVector &v, double alpha=1.);

/* *** Matrix-Vector operations ************************************************/
void CRSmatVec(DataVector &res, CRSMatrix &A, DataVector u);
void matVec(DataVector &res, Matrix &A, DataVector u);


#endif


