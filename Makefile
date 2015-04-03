CC=mpixlcxx
CFLAGS= -g -I./contrib
OMP=
LIBS=
cppfiles = $(shell ls src/*.cpp) $(shell ls contrib/*.cpp)
hfiles = $(shell ls src/*.h)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)
	$(CC) $(CFLAGS) $(ofiles) $(LIBS) -o coupledCellsModel

$(ofiles): %.o: %.cpp $(hfiles)
	$(CC) $(CFLAGS) $(OMP) $(LIBS) -o $@ -c $<

clean:
	rm -f src/*.o
	rm -f contrib/*.o
	rm -f coupledCellsModel

