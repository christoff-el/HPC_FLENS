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

all: FLENSDataVector.o cg_nompi_blas.o cg_mpi_blas.o gs_nompi_blas.o gs_mpi_blas.o FLENS_read_write_methods.o CRS_sort.o flens_supl.a

FLENSDataVector.o: FLENSDataVector.cpp FLENSDataVector.h
	$(CC) $(CPPFLAGS)  FLENSDataVector.cpp
flens2funk.o: FLENSDataVector.o flens2funk.cpp flens2funk.h
	$(CC) $(CPPFLAGS)  flens2funk.cpp
funk2flens.o: FLENSDataVector.o funk2flens.cpp funk2flens.h
	$(CC) $(CPPFLAGS)  funk2flens.cpp

cg_nompi_blas.o: cg_nompi_blas.cpp cg_nompi_blas.h
	$(CC) $(CPPFLAGS)  cg_nompi_blas.cpp
gs_nompi_blas.o: gs_nompi_blas.cpp gs_nompi_blas.h
	$(CC) $(CPPFLAGS)  gs_nompi_blas.cpp
cg_mpi_blas.o: cg_mpi_blas.cpp cg_mpi_blas.h
	$(CC) $(CPPFLAGS)  cg_mpi_blas.cpp
gs_mpi_blas.o: gs_mpi_blas.cpp gs_mpi_blas.h
	$(CC) $(CPPFLAGS)  gs_mpi_blas.cpp
wrappers.o: FLENSDataVector.o funk2flens.o flens2funk.o gs_mpi_blas.o gs_nompi_blas.o cg_mpi_blas.o cg_nompi_blas.o wrappers.cpp wrappers.h
	$(CC) $(CPPFLAGS)  wrappers.cpp
FLENS_read_write_methods.o: FLENS_read_write_methods.cpp FLENS_read_write_methods.h
	$(CC) $(CPPFLAGS)  FLENS_read_write_methods.cpp
CRS_sort.o: CRS_sort.cpp CRS_sort.h
	$(CC) $(CPPFLAGS)  CRS_sort.cpp 

flens_supl.a: FLENSDataVector.o cg_nompi_blas.o cg_mpi_blas.o gs_nompi_blas.o gs_mpi_blas.o FLENS_read_write_methods.o CRS_sort.o
	ar rcu libFlens.a FLENSDataVector.o cg_nompi_blas.o cg_mpi_blas.o gs_nompi_blas.o gs_mpi_blas.o FLENS_read_write_methods.o CRS_sort.o

clean:
	rm *.o libFlens.a
