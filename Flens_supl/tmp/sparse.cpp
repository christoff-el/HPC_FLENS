#include <iostream>

#include <flens/flens.cxx>

using namespace flens;
using namespace std;

class FLvNonMPI{};
class FLvTypeI{};
class FLvTypeII{};

template <typename VTYPE = FLvNonMPI>
struct FLENSDataVector
		: public DenseVector<Array<double> >
{
	typedef double	ElementType;
    typedef int   	IndexType;

	explicit
	FLENSDataVector(int n)
		:	DenseVector<Array<double> >(n)
		{
		cout<<"hi"<<endl;}
		
	explicit
	FLENSDataVector(FLENSDataVector &rhs)
		: DenseVector<Array<double> >(rhs)
		{}
	
	//DenseVector<Array<double> > vals;
    int uhoh=1;
    
    	
    
    
    


};



namespace flens {namespace blas {
void
copyer(FLENSDataVector<FLvTypeI> &a, FLENSDataVector<FLvTypeI> &b){
cout<<"Hahaha"<<endl;

//DenseVector<Array<double> > *tmpa = &a;
//DenseVector<Array<double> > *tmpb = &b;

blas::copy(*static_cast<DenseVector<Array<double> > *>(&b), *static_cast<DenseVector<Array<double> > *>(&a));}
    }}

int main() {

	typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordRowColCmp, IndexBase>  Coord;
    
    FLENSDataVector<> a(5);
    a(1)=99;
    FLENSDataVector<FLvTypeI> b(5);
    b(2)=101;
    flens::blas::copyer(a,b);
    
    cout<<a<<endl;
    cout<<b<<endl;
    
    
    /*FLENSDataVector b(5);
    b(2)=99;
    b(3)=1;
    
	FLENSDataVector a(5);
	
	flens::blas::copy(a,b);
	
	flens::GeCoordMatrix<Coord>  A_coord(4,4);
	A_coord(3,3)+=21;
	flens::GeCoordMatrix<CoordStorage<double> > u(3,3);
	u(3,3) += 808;
	cout<<u<<endl;
	
	flens::GeCRSMatrix<flens::CRS<double, flens::IndexOptions<int, 1>> > A = A_coord;
	
	cout<<A<<endl;
	DenseVector<Array<double, IndexBase> > q(10);
	q(0)=20;
	A.engine().values() = q;
	cout<<A<<endl;
	q(0)=10101;
	cout<<A<<endl;
	
	//blas::copy(b,a);
	
	cout << a << endl;
	cout << a.uhoh << endl;
	
	DenseVector<Array<double> > p(5);
	p(5) = 101;
	cout << p.length() << " " << p(5) << endl;
	p.resize(10);
	cout << p(5) << " " << p(6) << endl;*/
	
	
	//cout << a.sayHi(1) << endl;

	return 0;
}