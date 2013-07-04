OS = $(shell uname -s)

LDLIBS  := -L../LinearAlgebra -lLinAlg -L../Fem -lFem -L../Flens_supl -lFlens

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

CPPFLAGS = -Wall -I. -I.. -I/opt/local/include/openmpi $(LIBS) -I/usr/local/FLENS -std=c++11
#CPPFLAGS = -Wall -I. -I.. $(LIBS)

.PHONY:   all clean

all: main-serial main-parallel 
main-serial: main-serial.cpp
	$(CC) $(LDLIBS) $(CPPFLAGS) main-serial.cpp -o main-serial
main-parallel: main-parallel.cpp
	$(CC) $(LDLIBS) $(CPPFLAGS) main-parallel.cpp -o main-parallel

clean:
	rm -f main-serial main-parallel *.o 
