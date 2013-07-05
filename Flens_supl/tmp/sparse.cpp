#include <iostream>

#include <flens/flens.cxx>

using namespace flens;
using namespace std;


struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	/*explicit
	FLENSDataVector(int n)
		: DenseVector<Array<double> >(n),
		  vals(n)
		{}*/
	
	DenseVector<Array<double> > vals;
    
    int sayHi(int n);


};

int
FLENSDataVector::sayHi(int n){return (*this)(n);}

int main() {

	typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;
    
    DenseVector<Array<double> > b(5);
    b(2)=99;
    
	FLENSDataVector a;
	
	blas::copy(b,a);
	
	cout << a << endl;
	
	//cout << a.sayHi(1) << endl;

	return 0;
}