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

CPPFLAGS = -Wall -O3 -I. -I.. -I/opt/local/include/openmpi -I/opt/local/include -I/usr/local/FLENS -std=c++11 -c

.PHONY: all clean

all: cg_nompi_blas.o flens2funk.o funk2flens.o flens_supl.a
flens2funk.o: flens2funk.cpp flens2funk.h
	$(CC) $(CPPFLAGS)  flens2funk.cpp
funk2flens.o: funk2flens.cpp funk2flens.h
	$(CC) $(CPPFLAGS)  funk2flens.cpp
cg_nompi_blas.o: flens2funk.o funk2flens.o cg_nompi_blas.cpp cg_nompi_blas.h
	$(CC) $(CPPFLAGS)  cg_nompi_blas.cpp
flens_supl.a: flens2funk.o funk2flens.o cg_nompi_blas.o 
	ar rcu libFlens.a flens2funk.o funk2flens.o cg_nompi_blas.o

clean:
	rm *.o libFlens.a
