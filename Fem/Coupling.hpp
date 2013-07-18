#ifndef COUPLING_H_
#define COUPLING_H_

#include <flens/flens.cxx>
#include <vector>


class Coupling{
public:

	typedef flens::DenseVector<flens::Array<int> > IVector;
	
	/* *** public Variables */
	IVector neighbourProcs;
	std::vector<IVector> boundaryNodes;
	std::vector<IVector> coupling2edges;
	IVector local2globalCrossPoints;
	IVector crossPointsBdryData;
	IVector crossPointsNumProcs;
	IVector colors;
	
	int maxColor, numCoupling, numCrossPoints;
	
	/* *** constructors */
	Coupling();

  	Coupling(IVector _neighbourProcs, std::vector<IVector> _boundaryNodes, 
  				std::vector<IVector> _coupling2edges, IVector _local2globalCrossPoints, 
  				IVector _crossPointsBdryData, IVector _crossPointsNumProcs, 
				IVector _colors, int _maxColor, int _numCoupling, int _numCrossPoints);

  	Coupling(const Coupling & rhs);
	
	/* *** destructor */
	~Coupling();
	
};

#endif

