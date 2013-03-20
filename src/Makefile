CC=mpixlcxx_r
DEBUG= -g
CFLAGS=-c $(DEBUG)
LDFLAGS=
SOURCES=main.cpp compute.cpp checkpoint.cpp rksuite.cpp update.cpp solver_rk.cpp debug_file.cpp topology.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=output/test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
run: 	
	mpirun -np 12 ./test
clean: 
		rm output/*.txt *.o $(EXECUTABLE) output/*.err output/*.out

