CC=mpixlcxx
CFLAGS= -pg 
OMP= -openmp
LIBS=
cppfiles = $(shell ls *.cpp)
hfiles = $(shell ls *.h)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)
	$(CC) $(CFLAGS) $(ofiles) $(LIBS) -o model

cvode: CFLAGS=-pg -DCVODE
cvode: LIBS= -I/bgp/local/pkg/sundials/2.5.0/include/ -L/bgp/local/pkg/sundials/2.5.0/lib -lsundials_cvode -lsundials_nvecserial
cvode: all

$(ofiles): %.o: %.cpp $(hfiles)
	$(CC) $(CFLAGS) $(OMP) $(LIBS) -o $@ -c $<

clean:
	rm -f $(ofiles) model *.txt gout.* analysis*.txt

run:
	mpirun -np 12 ./model
