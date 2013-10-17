P=main
OBJECTS=cupgamelib.o
# CFLAGS=-Wall -O2 -funroll-loops -mfpmath=sse -m64 -Ofast -msse4 -march=native -fstrict-aliasing -fopenmp
CFLAGS=-Wall -O3 -mfpmath=sse -m64 -Ofast -msse4 -march=native -fstrict-aliasing -fopenmp -funroll-loops
LDLIBS=-lm -lrt $(OPTLIBS)
# CC=clang
CC=c99

$(P): $(OBJECTS) 
