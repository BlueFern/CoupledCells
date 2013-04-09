CC=mpixlcxx
CFLAGS=
LIBS=-lmass
INC = -I /bgp/local/pkg/sundials/2.5.0/lib/
cppfiles = $(shell ls *.cpp)
hfiles = $(shell ls *.h)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)
	$(CC) $(ofiles) $(LIBS) -o model

cvode: CFLAGS=-DCVODE
cvode: LIBS=-lmass -lsundials_cvode -lsundials_nvecserial
cvode: all

mac: CC=mpicxx
mac: LIBS=
mac: all

cvmac: CC=mpicxx
cvmac: CFLAGS=-DCVODE
cvmac: LIBS=-lsundials_cvode -lsundials_nvecserial
cvmac: all

$(ofiles): %.o: %.cpp $(hfiles)
	$(CC) -o $@ -c $<

clean:
	rm -f $(ofiles) model

