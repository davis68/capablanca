all: main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	mpicxx -Os -o capablanca main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	rm -f *.o *.gch

optimized: main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	mpicxx -Os -o capablanca main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o

unoptimized: main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o
	mpicxx -o capablanca main.o globals.o readxyz.o neighbors.o process.o error.o input.o statistics.o

debug: maing.o globalsg.o readxyzg.o neighborsg.o processg.o errorg.o inputg.o statisticsg.o
	mpicxx -g -o capablanca maing.o globalsg.o readxyzg.o neighborsg.o processg.o errorg.o inputg.o statisticsg.o

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

globalsg.o: globals.cpp
	mpicxx -g -o globalsg.o -c globals.cpp

maing.o: main.cpp
	mpicxx -g -o maing.o -c main.cpp

readxyzg.o: readxyz.cpp
	mpicxx -g -o readxyzg.o -c readxyz.cpp

neighborsg.o: neighbors.cpp
	mpicxx -g -o neighborsg.o -c neighbors.cpp

processg.o: process.cpp
	mpicxx -g -o processg.o -I /usr/include/boost -c process.cpp

errorg.o: error.cpp
	mpicxx -g -o errorg.o -c error.cpp

inputg.o: input.cpp
	mpicxx -g -o inputg.o -c input.cpp

statisticsg.o: statistics.cpp
	mpicxx -g -o statisticsg.o -c statistics.cpp

clean:
	rm -f *.o *.gch capablanca core* *.xyz_temp *.dat_temp surface*.xyz N=*

