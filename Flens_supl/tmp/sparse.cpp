#include <iostream>

#include <flens/flens.cxx>

using namespace flens;
using namespace std;


struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	explicit
	FLENSDataVector(int n)
		:	DenseVector<Array<double> >(n)
		{}
		
	explicit
	FLENSDataVector(FLENSDataVector &rhs)
		: DenseVector<Array<double> >(rhs)
		{}
	
	//DenseVector<Array<double> > vals;
    
    int sayHi(int n);


};

int
FLENSDataVector::sayHi(int n){return (*this)(n);}

int main() {

	typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;
    
    FLENSDataVector b(5);
    b(2)=99;
    
	FLENSDataVector a(b);
	
	//blas::copy(b,a);
	
	cout << a << endl;
	
	//cout << a.sayHi(1) << endl;

	return 0;
}