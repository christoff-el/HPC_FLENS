OS = $(shell uname -s)

.PHONY: all clean

ifeq ($(OS),Darwin)
  CC  := mpic++
  #CC  := openmpic++
endif

ifeq ($(OS),SunOS)
  CC  := mpiCC
endif

ifeq ($(OS),Linux)
  CC  := mpicxx
endif

CPPFLAGS = -Wall -O3 -I. -I.. -I/opt/local/include/openmpi -I/opt/local/include -std=c++11 -c

.PHONY: all clean

all: Vector.o Matrix.o CRSMatrix.o libLinAlg.a
Vector.o: Vector.cpp Vector.hpp 
	$(CC) $(CPPFLAGS)  Vector.cpp
Matrix.o: Vector.o Matrix.cpp Matrix.hpp 
	$(CC) $(CPPFLAGS)  Matrix.cpp
CRSMatrix.o: Vector.o Matrix.o CRSMatrix.cpp CRSMatrix.hpp 
	$(CC) $(CPPFLAGS)  CRSMatrix.cpp
libLinAlg.a: Vector.o Matrix.o CRSMatrix.o
	ar rcu libLinAlg.a Vector.o Matrix.o CRSMatrix.o

clean:
	rm *.o libLinAlg.a
