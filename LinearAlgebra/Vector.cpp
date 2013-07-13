#include <fstream>
#include <iostream>
#include "Vector.hpp"



/* ******************************************************************* */
/* ************************    IndexVector     *********************** */
/* ******************************************************************* */

/* *** Constructor & Destructor *************************************/
IndexVector::IndexVector(int length){
	assert(length>=0);
	_length=length;
	if(_length>0) {
		_data=new int[_length];
		memset(_data,0,_length*sizeof(int));
	}else{
		_data=0;
	}	
}

IndexVector::IndexVector(){
	_length=0;
	_data=0;
}

IndexVector::IndexVector(const IndexVector &rhs){
		_length=rhs.length();
		if(_length>0) {
			_data = new int[_length];
			memcpy(_data,rhs.data(),_length*sizeof(int) );
		}	else{
			_data=0;
		}
}

IndexVector::~IndexVector(){
	if(_data!=0) delete [] _data;
	_data=0;
}


/* *** getter and setter mehtods ***********************************/
int IndexVector::length() const{
	return _length;
}

void IndexVector::set(int index,int num,const int* val)
{	
	assert(index+num<=_length && index>=0);
	memcpy(_data+index,val, num*sizeof(int) );
}

void IndexVector::init(int val){	
	memset(_data,val,_length*sizeof(int));
}

int* IndexVector::data() const{
	return _data;
}

/* *** operators  *************************************************/
IndexVector & IndexVector::operator=(const IndexVector &rhs){
	if(_length!=rhs.length() )
		resize(rhs.length());
	
	memcpy( _data, rhs.data(), _length*sizeof(int) );
	return *this;
}

int &  IndexVector::operator()(const int n) const{
	assert(n<_length && n>=0);
    return _data[n];
}

/* *** math operations  *******************************************/
int IndexVector::dot(IndexVector &x){
	assert(_length==x.length());
	int tmp=0;
	for(int k=0;k<_length;k++) tmp+= _data[k]*x(k);
	return tmp;
}

int IndexVector::max(){
	int max=0;
	for(int k=0;k<_length;k++)
		if(_data[k]>max) max=_data[k];
		
	return max;
	
}

void IndexVector::add(IndexVector &x, int alpha){
	assert(_length==x.length());
	for(int k=0;k<_length;k++) _data[k]+=alpha*x(k);
}

void IndexVector::mul(int alpha){
	for(int k=0;k<_length;k++) _data[k]*=alpha;
}

/* *** write and read *********************************************/
void IndexVector::write(std::string filename){
	int j;
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open()){
		for(j=0;j<_length;j++){
			f << _data[j] << std::endl;
		}	
		f.close();
	}
}

void IndexVector::read(std::string filename){
	int j;
	FILE * in;
	if( (in = fopen(filename.c_str(),"r")) ==NULL){
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d",&j)==0){
		_length=0;
		_data=0;
	}else{
		resize(j);
		for(j=0;j<_length;j++){
			fscanf(in,"%d",&_data[j]);
		}
	}	
}
/* *** other functions ********************************************/
void IndexVector::resize(int length){
	if(_length!=length){
		_length=length;
		if(_data!=0) delete []_data;
		if(_length>0){
			_data=new int[_length];
			memset(_data,0,_length*sizeof(int));
		}else{
			_data=0;
		}
	}else{
		memset(_data,0,_length*sizeof(int));	
	}
}


/* ******************************************************************* */
/* ************************       Vector       *********************** */
/* ******************************************************************* */

/* *** Constructor & Destructor *************************************/
Vector::Vector(int length){
	assert(length>=0);
	_length=length;
	if(_length>0) {
		_data=new double[_length];
		memset(_data,0,_length*sizeof(double));
	}else{
		_data=0;
	}
}

Vector::Vector(){
	_length=0;
	_data=0;
}

Vector::Vector(const Vector &rhs){
	_length=rhs.length();
	if(_length>0) {
		_data = new double[_length];
		memcpy(_data,rhs.data(),_length*sizeof(double) );
	}	else{
		_data=0;
	}
}

Vector::~Vector(){
	if(_data!=0) delete [] _data;
	_data=0;
}

/* *** getter and setter mehtods ***********************************/
int Vector::length() const{
	return _length;
}

void Vector::set(int index,int num,double * val){	
	assert(index+num<=_length && index>=0);
	memcpy(_data+index,val, num*sizeof(double) );
}

void Vector::init(double val){	
	memset(_data,val,_length*sizeof(double));
}


double* Vector::data() const{
	return _data;
}
/* *** operators  *************************************************/
Vector& Vector::operator=(const Vector &rhs){
	if(_length!=rhs.length() )
		resize(rhs.length());
	
	memcpy( _data, rhs.data(), _length*sizeof(double));
	return *this;
}

double &  Vector::operator()(const int n) const {
	assert(n<_length && n>=0);
    return _data[n];
}

/* *** math operations  *******************************************/
double Vector::dot(Vector &x){
	assert(_length==x.length());
	double tmp=0;
	for(int k=0;k<_length;k++) tmp+= _data[k]*x(k);
	return tmp;
}

double Vector::max(){
	double max=0;
	for(int k=0;k<_length;k++)
		if(_data[k]>max) max=_data[k];
		
	return max;
}

void Vector::add(Vector &x, double alpha){
	assert(_length==x.length());
	for(int k=0;k<_length;k++) _data[k]+= alpha*x(k);
}

void Vector::mul(double alpha){
	for(int k=0;k<_length;k++) _data[k]*=alpha;
}

/* *** write and read *********************************************/
void Vector::write(std::string filename){
	int j;
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open()){
		for(j=0;j<_length;j++){
			f << _data[j] << std::endl;
		}
		f.close();
	}
}

void Vector::read(std::string filename){
	int j;
	FILE * in;
	if( (in = fopen(filename.c_str(),"r")) ==NULL){
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d",&j)==0){
		_length=0;
		_data=0;
	}else{
		resize(j);
		for(j=0;j<_length;j++){
			fscanf(in,"%lf",&_data[j]);
		}
	}	
}
/* *** other functions ********************************************/
void Vector::resize(int length){
	if(_length!=length){
		_length=length;
		if(_data!=0) delete [] _data;
		if(_length>0){
			_data=new double[_length];
			memset(_data,0,_length*sizeof(double));
		}else{
			_data=0;
		}
	} else{
		memset(_data,0,_length*sizeof(double));	
	}
}
