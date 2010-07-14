all: main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	mpicxx -o capablanca main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	rm -f *.o *.gch

optimized: main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	mpicxx -Os -o capablanca main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o

globals.o: globals.cpp
	mpicxx -c globals.cpp

main.o: main.cpp
	mpicxx -c main.cpp

readxyz.o: readxyz.cpp
	mpicxx -c readxyz.cpp

neighbors.o: neighbors.cpp
	mpicxx -c neighbors.cpp

process.o: process.cpp
	mpicxx -I /usr/include/boost -c process.cpp

error.o: error.cpp
	mpicxx -c error.cpp

input.o: input.cpp
	mpicxx -c input.cpp

statistics.o: statistics.cpp
	mpicxx -c statistics.cpp

clean:
	rm -f *.o *.gch capablanca core* *.xyz_temp *.dat_temp surface*.xyz N=*

