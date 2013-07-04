#ifndef MATRIX_H_
#define MATRIX_H_

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string.h>
#include <memory>

#include <LinearAlgebra/Vector.hpp>
//#include <LinearAlgebra/DataVector.hpp>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

class IndexMatrix{
	public:
		/* *** Constructor & Destructor *********************************************************/
		IndexMatrix();
		IndexMatrix(int rows, int cols);
		IndexMatrix(const IndexMatrix &rhs);
		~IndexMatrix();
		/* *** getter and setter functions ******************************************************/
		int numRows() const;
		int numCols() const;
		int* data() const;
		void set(int row, int numRows, int col, int numCols, int* val);

		/* *** operators *************************************************************************/
		IndexMatrix & operator =(const IndexMatrix& rhs);
		int &  operator()(const int m, const int n) const;
		IndexVector & matVec(IndexVector &res, IndexVector &x);
		
		/* *** write and read ********************************************************************/
		void write(std::string filename);
		void read(std::string filename);
		/* *** other functions *******************************************************************/
		void resize(int numcols,int numrows);
	
	/* *** private variables********************************************************************/
	private:
		int _numRows;
		int _numCols;
		int* _data;	
};



class Matrix{
public:
		/* *** constructor and destructor *******************************************************/
		Matrix();
		Matrix(int rows, int cols);
		Matrix(const Matrix &rhs);
		~Matrix();
		/* *** getter and setter functions *****************************************************/
		int numRows() const;
		int numCols() const;
		double* data() const;
		void set(int row, int numRows, int col, int numCols, double* val);
		/* *** operators ***********************************************************************/
		Matrix & operator =(const Matrix& rhs);
		double &  operator()(const int m, const int n) const;
		Vector & matVec(Vector &res, Vector &x);
		/* *** write and read ******************************************************************/
		void write(std::string filename);
		void read(std::string filename);
		/* *** other functions *****************************************************************/
		void resize(int numcols,int numrows);
		
	/* *** private variables *****************************************************************/
	private:
		int _numRows;
		int _numCols;
		double* _data;
};

#endif
