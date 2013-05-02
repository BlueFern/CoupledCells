CC=mpixlcxx
CFLAGS= 
OMP= -openmp
LIBS=
cppfiles = $(shell ls *.cpp)
hfiles = $(shell ls *.h)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)
	$(CC) $(CFLAGS) $(ofiles) $(LIBS) -o model

cvode: CFLAGS=-pg -DCVODE
cvode: LIBS= -I /intrepid-fs0/users/ccells/persistent/sundials/include -L /intrepid-fs0/users/ccells/persistent/sundials/lib -lsundials_cvode -lsundials_nvecserial
cvode: all

$(ofiles): %.o: %.cpp $(hfiles)
	$(CC) $(CFLAGS) $(OMP) $(LIBS) -o $@ -c $<

clean:
	rm -f $(ofiles) model *.txt gout.* analysis*.txt
