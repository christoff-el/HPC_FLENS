#include <iostream>
#include <string.h>
#include "FLENS_read_write_methods.h"
#include "CRS_sort.h"


using namespace std;

int
main()
{	
	flens::GeMatrix<flens::FullStorage<int> > A(4,4), B, C;
	flens::DenseVector<flens::Array<int> > rows(6), cols(6);
	flens::DenseVector<flens::Array<double> > V(5), vals(6), v;
	vals = 1;
	cout << "vals= " << vals << endl;
	V = 33;
	A = 6;
	//double &a = A(2,2);
	//a = 333;
	write(A, "test.txt");
	write(V, "test_vec.txt");
	read(C, "../examples/input/Square_four_domains/coordinates.dat");
	read(B, "test.txt");

	cout << C << endl << endl << endl;
	cout << C.numRows() << endl << C.numCols();

	return 0;
};