OS = $(shell uname -s)

LDFLAGS  := -L../Fem -L../Flens_supl
LDLIBS	 := -lFem -lFlens

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

CPPFLAGS = -Wall -I. -I.. -I../Flens_supl -I../Fem -I/opt/local/include/openmpi -I/usr/local/FLENS -std=c++11

.PHONY:   all clean

all: main-serial main-parallel 
main-serial: main-serial.cpp
	$(CC) $(LDFLAGS) $(CPPFLAGS) main-serial.cpp -o main-serial $(LDLIBS)
main-parallel: main-parallel.cpp
	$(CC) $(LDFLAGS) $(CPPFLAGS) main-parallel.cpp -o main-parallel $(LDLIBS)

clean:
	rm -f main-serial main-parallel *.o 
