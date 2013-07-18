OS = $(shell uname -s)
PC = $(shell uname -n)

LDFLAGS  := -L../Fem -L../Flens_supl -L../LinearAlgebra
LDLIBS	 := -lFem -lFlens -lLinAlg

ifeq ($(PC),pacioli)
  IFLENS   :=-I../../FLENS
else
  IFLENS   := -I/usr/local/FLENS
endif

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

CPPFLAGS = -Wall -I. -I.. -I../Flens_supl -I../Fem -I/opt/local/include/openmpi $(IFLENS) -std=c++11

.PHONY:   all clean

all: main-serial main-parallel 
main-serial: main-serial.cpp
	$(CC) $(LDFLAGS) $(CPPFLAGS) main-serial.cpp -o main-serial $(LDLIBS)
main-parallel: main-parallel.cpp
	$(CC) $(LDFLAGS) $(CPPFLAGS) main-parallel.cpp -o main-parallel $(LDLIBS)

clean:
	rm -f main-serial main-parallel *.o 
