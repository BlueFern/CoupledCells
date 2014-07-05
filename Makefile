CC=mpixlcxx
CFLAGS= -g
OMP= 
LIBS=
cppfiles = $(shell ls src/*.cpp)
hfiles = $(shell ls src/*.h)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)
	$(CC) $(CFLAGS) $(ofiles) $(LIBS) -o model

cvode: CFLAGS= -DCVODE
cvode: LIBS= -I/bgp/local/pkg/sundials/2.5.0/include/ -L/bgp/local/pkg/sundials/2.5.0/lib -lsundials_cvode -lsundials_nvecserial
cvode: all

$(ofiles): %.o: %.cpp $(hfiles)
	$(CC) $(CFLAGS) $(OMP) $(LIBS) -o $@ -c $<

clean:
	rm -f $(ofiles) model *.txt gout.* analysis*.txt

run:
	mpirun -np 12 ./model
