all: main.o globals.o readxyz.o nn.o process.o error.o input.o statistics.o
	mpicxx -o capablanca main.o globals.o readxyz.o nn.o process.o error.o input.o statistics.o

optimized:
	mpicxx -Os -o capablanca main.o globals.o readxyz.o nn.o process.o error.o input.o statistics.o

globals.o: globals.cc
	mpicxx -c globals.cc

main.o: main.cpp
	mpicxx -c main.cpp

readxyz.o: readxyz.cc
	mpicxx -c readxyz.cc

nn.o: findNnew.cc
	mpicxx -o nn.o -c findNnew.cc

process.o: process.cc
	mpicxx -c process.cc

error.o: error.cpp
	mpicxx -c error.cpp

input.o: input.cpp
	mpicxx -c input.cpp

statistics.o: statistics.cc
	mpicxx -c statistics.cc

clean:
	rm -f *.o *.gch capablanca core* *.xyz_temp *.dat_temp

