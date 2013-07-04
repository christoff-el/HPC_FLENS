#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <map>
#include "Matrix.hpp"



/* ***************************************************************************************** */
/* ********************************    IndexMatrix     ************************************* */
/* ***************************************************************************************** */

/* ***Constructor & Destructor ***************************************************************/
IndexMatrix::IndexMatrix(int rows, int cols){
	_numRows=rows;
	_numCols=cols;
	if(_numRows*_numCols>0){ 
		_data = new int[_numRows*_numCols];
		memset(_data,0,_numRows*_numCols*sizeof(int));
	}else{
		_data=0;
	}
}

IndexMatrix::IndexMatrix(){
	_numRows=0;
	_numCols=0;
	_data=0;
}

IndexMatrix::IndexMatrix(const IndexMatrix &rhs){
		_numRows=rhs.numRows();
		_numCols=rhs.numCols();
		if(_numRows>0 && _numCols>0) {
			_data = new int[_numRows*_numCols];
			memcpy(_data,rhs.data(),_numRows*_numCols*sizeof(int) );
		}else{
			_data=0;
		}
}

IndexMatrix::~IndexMatrix(){
	if(_data!=0) delete []_data;
	_data=0;
}


/* *** getter and setter methods ***************************************************************/
int IndexMatrix::numRows() const{
	return _numRows;
}

int IndexMatrix::numCols() const{
	return _numCols;
}

int* IndexMatrix::data() const{
	return _data;
}

void IndexMatrix::set(int row, int numRows, int col, int numCols, int* val){
	assert(row+numRows <=_numRows && col+numCols <= _numCols);
	for(int k=0;k<numCols;k++){
		memcpy(_data+(k+col)*_numRows+row, (val+k*numRows), numRows*sizeof(int) );
	}
}

/* *** operators *******************************************************************************/
IndexMatrix & IndexMatrix::operator=(const IndexMatrix &rhs){
	if(_numRows!=rhs.numRows() || _numCols!=rhs.numCols())
		resize(rhs.numRows(),rhs.numCols());
	
	memcpy(_data, rhs.data(), _numRows*_numCols*sizeof(int) );
	return *this;
}

int & IndexMatrix::operator()(const int m,const int n) const{
	assert(m<_numRows && m>=0 && n<_numCols && n>=0);
    return _data[m+n*_numRows];
}

IndexVector& IndexMatrix::matVec(IndexVector &res, IndexVector &x){
	assert( _numCols == x.length() && _numRows==res.length() );

	for(int j=0;j<_numRows;j++){
		for(int k=0;k<_numCols;k++){
			res(j) += _data[k*_numRows+j]*x(k);
		}
	}
	return res;
}

/* *** write and read methods ******************************************************************/
void IndexMatrix::write(std::string filename)
{
	int j,k;
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open()){
		for(k=0;k<_numRows;k++){
			for(j=0;j<_numCols;j++){
				f << _data[_numRows*j+k] << "  ";
			}
			f << std::endl;
		}
		f << std::endl;
		f.close();
	}	
}

void IndexMatrix::read(std::string filename)
{
	int j,k;
	FILE * in;
	if( (in = fopen(filename.c_str(),"r")) ==NULL){
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d %d",&j, &k)==0){
		_numRows=0;
		_numCols=0;
		_data=0;
	}else{
		resize(j,k);
		for(int l=0;l<_numRows;l++)
			for(int m=0;m<_numCols;m++){
				fscanf(in,"%d",&_data[_numRows*m+l]);
			}
	}	
}
/* *** other methods ************************************************************************/
void IndexMatrix::resize(int numrows,int numcols)
{
	if(_numRows!=numrows || _numCols!=numcols){
		_numRows = numrows;
		_numCols = numcols;
		if(_data!=0) delete []_data;
		if(_numRows*_numCols>0){
			_data=new int[_numRows*_numCols];
			memset(_data,0,_numRows*_numCols*sizeof(int));
		}else{
			_data=0;
		}
	}
}


/* *************************************************************************************** */
/* *********************************    Matrix     *************************************** */
/* *************************************************************************************** */

/* *** Constructors & Destructor ***********************************************************/
Matrix::Matrix(int rows, int cols){
	_numRows=rows;
	_numCols=cols;
	if(_numRows*_numCols>0){ 
		_data = new double[_numRows*_numCols];
		memset(_data,0,_numRows*_numCols*sizeof(double));
	}else{
		_data=0;
	}
}

Matrix::Matrix(){
	_numRows=0;
	_numCols=0;
	_data=0;
}

Matrix::Matrix(const Matrix &rhs){
	
	_numRows=rhs.numRows();
	_numCols=rhs.numCols();
	if(_numRows>0 && _numCols>0) {
		_data = new double[_numRows*_numCols];
		memcpy(_data,rhs.data(),_numRows*_numCols*sizeof(double) );
	}else{
		_data=0;
	}
}

Matrix::~Matrix(){
	if(_data!=0) delete []_data;
	_data=0;
}

/* *** operators **************************************************************************/
Matrix & Matrix::operator=(const Matrix &rhs){
	if(_numRows!=rhs.numRows() || _numCols!=rhs.numCols())
		resize(rhs.numRows(),rhs.numCols());
	
	memcpy(_data, rhs.data(), _numRows*_numCols*sizeof(double) );
	return *this;
}

double &  Matrix::operator()(const int m,const int n) const {
	assert(m<_numRows && m>=0 && n<_numCols && n>=0);
    return _data[m+n*_numRows];
}

Vector & Matrix::matVec(Vector &res, Vector &x){
	assert( _numCols==x.length() && _numRows==res.length() );

	for(int j=0;j<_numRows;j++){
		for(int k=0;k<_numCols;k++){
			res(j) += _data[k*_numRows+j]*x(k);
		}
	}
	return res;
}

/* *** getter and setter methods ***********************************************************/
int Matrix::numRows() const {
	return _numRows;
}

int Matrix::numCols() const{
	return _numCols;
}


double* Matrix::data() const{
	return _data;
}

void Matrix::set(int row, int numRows, int col, int numCols, double* val){
	assert(row+numRows <=_numRows && col+numCols <= _numCols);
	for(int k=0;k<numCols;k++){
		memcpy(_data+(k+col)*_numRows+row, (val+k*numRows), numRows*sizeof(double) );
	}
}

/* *** write and read methods ***********************************************************/
void Matrix::write(std::string filename)
{
	int j,k;
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open()){
		for(k=0;k<_numRows;k++){
			for(j=0;j<_numCols;j++){
				f.width(4);
				f << _data[_numRows*j+k] << "  ";
			}
			f << std::endl;
		}
		f << std::endl;
		f.close();
	}	
}

void Matrix::read(std::string filename)
{
	int j,k;
	FILE * in;
	if( (in = fopen(filename.c_str(),"r")) ==NULL){
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d %d",&j, &k)==0){
		_numRows=0;
		_numCols=0;
		_data=0;
	}else{
		resize(j,k);
		for(int l=0;l<_numRows;l++)
			for(int m=0;m<_numCols;m++){
				fscanf(in,"%lf",&_data[_numRows*m+l]);
			}
	}
}

/* *** other methods ***********************************************************/
void Matrix::resize(int numrows,int numcols){
	
	if(_numRows!=numrows || _numCols!=numcols){
		_numRows = numrows;
		_numCols = numcols;
		if(_data!=0) delete []_data;
		if(_numRows*_numCols>0){
			_data=new double[_numRows*_numCols];
			memset(_data,0,_numRows*_numCols*sizeof(double));
		}else{
			_data=0;
		}
	}
}

