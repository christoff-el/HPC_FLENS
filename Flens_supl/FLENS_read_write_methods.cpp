#ifndef FLENS_READ_WRITE_METHODS_CPP
#define FLENS_READ_WRITE_METHODS_CPP 1

#include "FLENS_read_write_methods.h"

namespace flens{

// Read Matrix
template <typename ElementType>
void
read(flens::GeMatrix<flens::FullStorage<ElementType> > &GeMatrix, std::string filename)
{
	typedef int IndexType;
	IndexType j,k;
	FILE * in;

	// Format check
	std::string format;
	if ((typeid(ElementType)==typeid(int))
		|| (typeid(ElementType)==typeid(unsigned int))
		|| (typeid(ElementType)==typeid(long int))
		|| (typeid(ElementType)==typeid(long unsigned int))) format = "%d";
	else format = "%lf";


	if( (in = fopen(filename.c_str(),"r")) ==NULL)
	{
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d %d",&j, &k)==0)
	{
		GeMatrix.resize(0,0);
	}
	else
	{ 
		GeMatrix.resize(j,k);
		for(IndexType r=1; r<=GeMatrix.numRows(); ++r)
		{
			for(IndexType c=1; c<=GeMatrix.numCols(); ++c)
			{
				fscanf(in, format.c_str(), &GeMatrix(r,c));
			}
		}
	}
};


// Read Vector
template <typename ElementType>
void
read(flens::DenseVector<flens::Array<ElementType> > &DVector, std::string filename)
{
	typedef int IndexType;
	IndexType j;

	// Format check
	std::string format;
	if ((typeid(ElementType)==typeid(int))
		|| (typeid(ElementType)==typeid(unsigned int))
		|| (typeid(ElementType)==typeid(long int))
		|| (typeid(ElementType)==typeid(long unsigned int))) format = "%d";
	else format = "%lf";


	FILE * in;
	if( (in = fopen(filename.c_str(),"r")) ==NULL)
	{
		perror("read: ");
		return;
	}
	if (fscanf(in,"%d",&j)==0)
	{
		DVector.resize(0);
	}
	else
	{
		DVector.resize(j);
		for(IndexType l=1; l<=DVector.length(); ++l)
		{
			fscanf(in, format.c_str(), &DVector(l));
		}
	}
};


// Write Matrix
template <typename ElementType>
void
write(flens::GeMatrix<flens::FullStorage<ElementType> > &GeMatrix, std::string filename)
{
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open())
	{
		f << GeMatrix;
		f.close();
	}
};

// Write Vector
template <typename ElementType>
void
write(flens::DenseVector<flens::Array<ElementType> > &DVector, std::string filename)
{
	typedef int IndexType;
	std::fstream f;
	
	f.open(filename.c_str(), std::ios::out);
	if(f.is_open())
	{
		for (IndexType i=1; i<=DVector.length(); ++i)
		{
			f << DVector(i) << std::endl;
		}
		f.close();
	}
};

}; //namespace flens

#endif	// FLENS_READ_WRITE_METHODS_CPP
