#include <algorithm>
#include <flens/flens.cxx>
#include <iostream>
#include "gs_nompi_blas.h"

using namespace flens;
using namespace std;

int
main()
{
	typedef GeMatrix<FullStorage<double, ColMajor> > GEMatrix;
	typedef DenseVector<Array<double>> DenseVector;

	GEMatrix A(3,3);
	//fillRandom(A);
	A = 135, 115, 5,
		23, 200, 33,
		1, 8, 18;
	DenseVector x(3), b(3), c(3);
	//fillRandom(x);
	x = 1, 2, 4;
	b = A*x;
	cout << "x=" << x << endl;
	x = 0;

	int It = gs_dense_nompi_blas(A, b, x, c, 30);

	cout << "It=" << It << endl; 
	cout << "A= " << A << endl;
	cout << "x= " << x << endl;
	cout << "b=" << b << endl;
	return 0;
};