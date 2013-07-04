#include "Coupling.hpp"

Coupling::Coupling() : neighbourProcs(), boundaryNodes(), coupling2edges(), local2globalCrossPoints(), crossPointsBdryData(), 
                       crossPointsNumProcs(), colors(), maxColor(0), numCoupling(0), numCrossPoints(0)
{	

}

Coupling::Coupling( IndexVector _neighbourProcs, std::vector<IndexVector> _boundaryNodes, std::vector<IndexVector> _coupling2edges,
          					IndexVector _local2globalCrossPoints, IndexVector _crossPointsBdryData, IndexVector _crossPointsNumProcs, 
										IndexVector _colors, int _maxColor, int _numCoupling, int _numCrossPoints):
										neighbourProcs(_neighbourProcs), boundaryNodes(_boundaryNodes), coupling2edges(_coupling2edges), 
										local2globalCrossPoints(_local2globalCrossPoints), crossPointsBdryData(_crossPointsBdryData),
										crossPointsNumProcs(_crossPointsNumProcs), colors(_colors), maxColor(_maxColor), numCoupling(_numCoupling),
										numCrossPoints(_numCrossPoints) 
{

}

Coupling::Coupling(const Coupling &rhs)
{
	neighbourProcs     		= rhs.neighbourProcs;
	boundaryNodes 			= rhs.boundaryNodes;
	coupling2edges 			= rhs.coupling2edges;
	local2globalCrossPoints	= rhs.local2globalCrossPoints;
	crossPointsBdryData 	= rhs.crossPointsBdryData;
	crossPointsNumProcs 	= rhs.crossPointsNumProcs;
	colors 					= rhs.colors;
	
	maxColor 				= rhs.maxColor;
	numCoupling 			= rhs.numCoupling;
	numCrossPoints 			= rhs.numCrossPoints;
}

Coupling::~Coupling(){
	
}
