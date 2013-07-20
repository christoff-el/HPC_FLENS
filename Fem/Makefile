OS = $(shell uname -s)
PC = $(shell uname -n)

LDFLAGS  := -L../Flens_supl
LDLIBS	 := -lFlens

ifeq ($(PC),pacioli)
  IFLENS   :=-I../../FLENS
else
  IFLENS   := -I/usr/local/FLENS
endif

.PHONY: all clean goto

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

CPPFLAGS = -Wall -O3 -I. -I.. -I/opt/local/include/openmpi -I/opt/local/include $(IFLENS) -std=c++11 -c


all: Coupling.o Mesh.o Fem.o libFem.a

goto: LDFLAGS  +=-L../../gotoins/lib
goto: LDLIBS   +=-lopenblas
goto: IFLENS   :=-I../../FLENS -DWITH_GOTOBLAS -I../../gotoins/include
goto: Coupling.o Mesh.o Fem.o libFem.a

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
Fem.o: Mesh.o Fem.cpp Fem.hpp
	$(CC) $(CPPFLAGS)  Fem.cpp


libFem.a:Coupling.o Mesh.o Fem.o
	ar rcu libFem.a  Coupling.o Mesh.o Fem.o


clean:
	rm *.o libFem.a 

