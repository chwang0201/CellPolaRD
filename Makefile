include ../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

CellPolaRD: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: CellPolaRD

run: all
	mpirun -np 4 ./gray_scott

.PHONY: clean all run

clean:
	rm -f *.o *~ core CellPolaRD

