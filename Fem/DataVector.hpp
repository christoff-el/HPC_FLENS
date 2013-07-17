#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_
#include <iostream>
#include <mpi.h>
#include <string>

#include "../LinearAlgebra/LinAlgHeader.hpp"
#include "Coupling.hpp"


enum vectorType { typeI, typeII, nonMPI};

class DataVector {
public:
    /* *** constructors ********************************************/
    DataVector(const Coupling &_coupling, const Vector &_values, 
               const vectorType _type=nonMPI);

    DataVector(const Coupling &_coupling,
               const int size=1 , const vectorType _type=nonMPI);

    DataVector(const DataVector &rhs);

    /* *** overloading operators ***********************************/
    DataVector& operator =(const DataVector &rhs);

    /* *** functions for conversion ********************************/
    void typeII_2_typeI();
    void typeI_2_typeII();

    /* *** functions for communication *****************************/
    void communicationCrossPoints();
    void communicationBoundaryNodes();

    /* *** other functions ****************************************/
    void resize(int N);

    /* *** producing output ***************************************/
    void
    writeData(int proc=0,std::string filename = "output/x");

    /* *** variables***********************************************/
    Vector values;
    vectorType type;
    const Coupling &coupling;
};



#endif

