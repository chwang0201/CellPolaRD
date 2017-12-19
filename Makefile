include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

CellPolaRD_Cha: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: CellPolaRD_Cha

run: all
	mpirun -np 4 ./CellPolaRD_Cha

.PHONY: clean all run

clean:
	rm -f *.o *~ core CellPolaRD_Cha

