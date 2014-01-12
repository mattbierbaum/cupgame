EXE=cupgame-nowait-mapgen
OBJECTS=cupgamelib.o
CFLAGS=-Wall -O3 -m64 -Ofast -march=native -fopenmp -D_POSIX_C_SOURCE=199309L
LDLIBS=-lm -lrt
CC=c99

$(EXE): $(OBJECTS) 

.PHONY: clean

clean:
	rm $(EXE) $(OBJECTS)
