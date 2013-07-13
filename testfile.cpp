#include <iostream>
#include <flens/flens.cxx>
using namespace flens;
using namespace std;

int
main()
{
	int data[4] = {1, 2, 3, 4};
	typedef DenseVector<Array<int> > Vek;
	Vek VV(4), A(4);
	const Vek &V = VV;
	typedef DenseVector<Array<int>::ConstView> Vec;
	const Underscore<int> _;
	Vec C = V(_(2,4));
	int* ptr = A.data()+1;
	int n = C(2);
	//A = B;
	data[1] = 333;
	cout << "V=" << V << endl;
	//cout << "B=" << B << endl;
	cout << "C=" << C << endl;	
	return 0;
}