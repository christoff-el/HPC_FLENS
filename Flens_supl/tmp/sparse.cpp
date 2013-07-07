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
    int uhoh=1;
    
    	
    
    
    


};



namespace flens {namespace blas {
void
copy(FLENSDataVector &a, FLENSDataVector &b){
cout<<"Hahaha"<<endl;

//DenseVector<Array<double> > *tmpa = &a;
//DenseVector<Array<double> > *tmpb = &b;

blas::copy(*static_cast<DenseVector<Array<double> > *>(&b), *static_cast<DenseVector<Array<double> > *>(&a));}
    }}

int main() {

	typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;
    
    FLENSDataVector b(5);
    b(2)=99;
    b(3)=1;
    
	FLENSDataVector a(5);
	
	flens::blas::copy(a,b);
	
	DenseVector<Array<double> > c(5);
	DenseVector<Array<double> > *d = &c;
	
	c(4) = 67;
	(*d)(5) = 76;
	
	cout << c << endl;
	cout << *d << endl;
	//blas::copy(b,a);
	
	cout << a << endl;
	cout << a.uhoh << endl;
	
	//cout << a.sayHi(1) << endl;

	return 0;
}