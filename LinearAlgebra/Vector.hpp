#ifndef VECTOR_H_
#define VECTOR_H_

#include <cassert>
#include <cstdlib>
#include <memory>
#include <string.h>
#include <string>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

class IndexVector{
	public:
		/* *** Constructor & Destructor *************************************/
		IndexVector();
		IndexVector(int length);
		IndexVector(const IndexVector &rhs);
		~IndexVector();
		/* *** getter and setter mehtods ***********************************/
		int length() const;
		int* data() const;
		void set(int index, int num, const int* val);
		void   init(int val);
		/* *** operators  *************************************************/
		IndexVector & operator=(const IndexVector &rhs);
		int & operator()(const int n) const;
		/* *** math operations  *******************************************/
		int dot(IndexVector &x);					// scalar product: returns this*x
		int max();									// maximum: returns max(this)
		void add(IndexVector &x, int alpha=1);		// this += alpha*x 
		void mul(int alpha);						// this *= alpha
		/* *** write and read *********************************************/
		void write(std::string filename);
		void read(std::string filename);
		/* *** other functions ********************************************/
		void resize(int length);
		
	/* *** private variables **********************************************/
	private:
		int _length;
		int* _data;
};


class Vector{
	
	public:
		/* *** Constructor & Destructor *************************************/
		Vector();
		Vector(int length);
		Vector(const Vector &rhs);
		~Vector();
		/* *** getter and setter mehtods ***********************************/
		int length() const;
		double* data() const;
		void   set(int index, int num, double* val);
		void   init(double val);
		/* *** operators  *************************************************/
		Vector & operator=(const Vector &rhs);
		double & operator()(const int n) const;
		/* *** math operations  *******************************************/
	  double dot(Vector &x);										// scalar product: returns this*x
		double max();															// maximum: returns max(this)
		void add(Vector &x, double alpha=1.);			// this += alpha*x 
		void mul(double alpha);										// this *= alpha
		/* *** write and read *********************************************/
		void write(std::string  filename);
		void read(std::string filename);
		/* *** other functions ********************************************/
		void resize(int length);
	/* *** private variables **********************************************/
	private:
		int _length;
		double* _data;	
};

#endif
