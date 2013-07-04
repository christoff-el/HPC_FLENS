#ifndef COUPLING_H_
#define COUPLING_H_

#include <vector>

#include "../LinearAlgebra/LinAlgHeader.hpp"


class Coupling{
public:
	/* *** public Variables */
	IndexVector neighbourProcs;
	std::vector<IndexVector> boundaryNodes;
	std::vector<IndexVector> coupling2edges;
	IndexVector local2globalCrossPoints;
	IndexVector crossPointsBdryData;
	IndexVector crossPointsNumProcs;
	IndexVector colors;
	
	int maxColor, numCoupling, numCrossPoints;
	
	/* *** constructors */
	Coupling();

  Coupling( IndexVector _neighbourProcs, std::vector<IndexVector> _boundaryNodes, std::vector<IndexVector> _coupling2edges,
            IndexVector _local2globalCrossPoints, IndexVector _crossPointsBdryData, IndexVector _crossPointsNumProcs, 
						IndexVector _colors, int _maxColor, int _numCoupling, int _numCrossPoints);

  Coupling(const Coupling & rhs);
	~Coupling();
	
};

#endif

