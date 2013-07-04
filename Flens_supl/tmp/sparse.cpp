#include <iostream>

#include <flens/flens.cxx>

using namespace flens;
using namespace std;


int main() {

	typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;
    
	GeCRSMatrix<CRS<double, IndexBase> >  B;
	
	//B(2,2)=4;
	
	cout << B << endl;


	return 0;
}