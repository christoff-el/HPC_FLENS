#ifndef CRSMATRIX_H_
#define CRSMATRIX_H_
#include <fstream>
#include <string>
#include <iostream>

#include "Matrix.hpp"
#include <string.h>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))


class CRSMatrix{
	
	public:
		/* *** Consturctor & Destructor */
		CRSMatrix();
		CRSMatrix(IndexVector &row, IndexVector &col,Vector &data);
	  CRSMatrix(const CRSMatrix &rhs);
	  ~CRSMatrix();

		/* *** getter & setter methods */
		int nonZeros() const;
	  int numCols() const;
	  int numRows() const;
		//double get(int row,int col) const;
		//double get(int k) const;
		double* data() const;
		int* colIndex() const;
		int* rowPtr() const;
	
		/* *** operators */
		CRSMatrix & operator=(const CRSMatrix &rhs);
		double  operator()(const int m, const int n);
		double &  operator()(const int m);
		Vector& matVec(Vector &res, Vector &x);
	
		/* *** convert routines */
		void CRS2full(Matrix &A);
		void full2CRS(Matrix &rhs);
	
		/* *** other routines */
		bool is_index(int row, int col);
		void writeCRS(std::string filename);
		void writeFull(std::string filename);
	
	
	private:
		int*     _rowptr;
		int*     _colindex;
		double*  _data;
		int      _cols;
		int      _rows;
		int      _nonZeros;
		void insort(int * index, double* data,int length);

};

class CRSIndexMatrix{
	
	public:
		/* *** Consturctor & Destructor */
		CRSIndexMatrix();
		CRSIndexMatrix(IndexVector &row, IndexVector &col,IndexVector &data);
	  CRSIndexMatrix(const CRSIndexMatrix &rhs);
	  ~CRSIndexMatrix();

		/* *** getter & setter methods */
		int nonZeros() const;
		//int get(int row,int col) const;
		//int get(int k) const;
	  int numCols() const;
	  int numRows() const;
		int* data() const;
		int* colIndex() const;
		int* rowPtr() const;
	
		/* *** operators */
		CRSIndexMatrix & operator=(const CRSIndexMatrix &rhs);
		int  operator()(const int m, const int n);
		int &  operator()(const int m);
		/* *** convert routines */
		void CRS2full(IndexMatrix &A);
		void full2CRS(IndexMatrix &rhs);
	
		/* *** other routines */
		bool is_index(int row, int col);
		void writeCRS(std::string filename);
		void writeFull(std::string filename);
	
	private:
		int*     _rowptr;
		int*     _colindex;
		int*  	 _data;
		int      _cols;
		int      _rows;
		int      _nonZeros;
		void insort(int * index, int* data,int length);

};

#endif	
