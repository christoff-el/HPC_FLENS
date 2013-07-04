#include "MathOperationsMPI.hpp"


/* *** Vector operations ************************************/
// compute global maximum
double max(DataVector& u)
{
	if(u.type==nonMPI) 
	  return u.values.max();
	
	if(u.type==typeII) 
	  u.typeII_2_typeI();
	
	double max_loc =  u.values.max();
	double max_gl;
	MPI::COMM_WORLD.Allreduce( &max_loc, &max_gl, 1,MPI::DOUBLE,MPI::MAX);
	
	return max_gl;
}

// compute global 2-norm 
double norm(DataVector& u)
{
	DataVector tmp(u);
	if(u.type==nonMPI) 
	  return sqrt( u.values.dot(tmp.values));
	
	if(tmp.type==typeI){ 
		tmp.typeI_2_typeII();
		return sqrt( dot(u,tmp) );	
	} else {
		tmp.typeII_2_typeI();
		return sqrt( dot(tmp, u) );	
	}
}

/* *** Vector-Vector operations *****************************/
// compute scalar product of typeI and typeII vector
double dot(DataVector &u, DataVector &v)
{
	if(u.type==nonMPI && v.type==nonMPI) 
	  return u.values.dot(v.values);
	
	// we only multiply typeI and typeII vectors
	assert(u.type != v.type);
	
	double value = u.values.dot(v.values);
	double buf=0;
	
	/* *** communication to add values from others procs */
  MPI::COMM_WORLD.Allreduce(&value, &buf, 1,MPI::DOUBLE,MPI::SUM);
  return buf;
}

// compute u += alpha*v
void add(DataVector &u, DataVector &v, double alpha)
{
	assert(u.type==v.type);
	u.values.add(v.values, alpha);
}

/* *** Matrix-Vector operations *****************************/
// CRS-Matrix * typeI-Vector
void CRSmatVec(DataVector &res, CRSMatrix &A, DataVector u)
{
	assert(u.type==nonMPI || u.type==typeI);
	if(u.type==nonMPI) 
	  res.type=nonMPI;
	else 
		res.type=typeII;
	
	A.matVec(res.values,u.values);
}

// Full-Matrix * typeI-Vector
void matVec(DataVector &res, Matrix &A, DataVector u)
{
	assert(u.type==nonMPI || u.type==typeI);
	if(u.type==nonMPI) 
	  res.type=nonMPI;
	else 
	  res.type=typeII;
	
	A.matVec(res.values,u.values);
}