#ifndef FLENS_READ_WRITE_METHODS_H
#define FLENS_READ_WRITE_METHODS_H 1

#include <fstream>
#include <string.h>
#include <typeinfo>
#include <flens/flens.cxx>

namespace flens{

// Read Matrix
template <typename ElementType>
void
read(flens::GeMatrix<flens::FullStorage<ElementType> > &GeMatrix, std::string filename);

// Read Vector
template <typename ElementType>
void
read(flens::DenseVector<flens::Array<ElementType> > &DVector, std::string filename);



// Write Matrix
template <typename ElementType>
void
write(flens::GeMatrix<flens::FullStorage<ElementType> > &GeMatrix, std::string filename);

// Write Vector
template <typename ElementType>
void
write(flens::DenseVector<flens::Array<ElementType> > &DVector, std::string filename);


}; // namespace flens

#include "FLENS_read_write_methods.cpp"
#endif // FLENS_READ_WRITE_METHODS_H