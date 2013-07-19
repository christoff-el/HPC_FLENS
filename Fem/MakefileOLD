OS = $(shell uname -s)

LDFLAGS  := -L../Flens_supl
LDLIBS	 := -lFlens

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

CPPFLAGS = -Wall -O3 -I. -I.. -I/opt/local/include/openmpi -I/opt/local/include -I/usr/local/FLENS -std=c++11 -c


all: Coupling.o DataVector.o MathOperationsMPI.o Mesh.o Solver.o Fem.o libFem.a
Coupling.o:  Coupling.cpp Coupling.hpp 
	$(CC) $(CPPFLAGS)  Coupling.cpp
DataVector.o:  Mesh.o DataVector.cpp DataVector.hpp 
	$(CC) $(CPPFLAGS)  DataVector.cpp
MathOperationsMPI.o: DataVector.o Mesh.o MathOperationsMPI.cpp MathOperationsMPI.hpp 
	$(CC) $(CPPFLAGS)  MathOperationsMPI.cpp
Mesh.o:  Coupling.o Mesh.cpp Mesh.hpp 
	$(CC) $(CPPFLAGS)  Mesh.cpp
Solver.o:  Mesh.o MathOperationsMPI.o DataVector.o Solver.cpp Solver.hpp 
	$(CC) $(CPPFLAGS)  Solver.cpp
Fem.o: Mesh.o DataVector.o Solver.o Fem.cpp Fem.hpp
	$(CC) $(CPPFLAGS) Fem.cpp


libFem.a:Coupling.o DataVector.o MathOperationsMPI.o Mesh.o Solver.o Fem.o
	ar rcu libFem.a  Coupling.o DataVector.o MathOperationsMPI.o Mesh.o Solver.o Fem.o


clean:
	rm *.o libFem.a 

