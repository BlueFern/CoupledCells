# MPICH Makefile for the MPI Tutorial examples
##### User configurable options #####

MPICH_HOME = /usr/local/

CC          = $(MPICH_HOME)/mpicc
CCC         = $(MPICH_HOME)/mpiCC
F77         = $(MPICH_HOME)/mpif77
F90         = $(MPICH_HOME)/mpif90
CLINKER     = $(MPICH_HOME)/mpicc
CCLINKER    = $(MPICH_HOME)/mpiCC
FLINKER     = $(MPICH_HOME)/mpif77
F90LINKER   = $(MPICH_HOME)/mpif90	  
OPTFLAGS    = -O3
MAKE        = make
SHELL       = /bin/sh
#

### End User configurable options ###

### to run, type "mpichrun -np 6 group_example" ###
### (assuming that you have already aliased mpichrun to $(MPICH_HOME)/mpirun)

CFLAGS  = $(OPTFLAGS)
FFLAGS =  $(OPTFLAGS)
LIBS =
FLIBS =
CLIBS =
EXECS = group_example comm_split_example

default: $(EXECS)

all: $(EXECS)

group_example: group_example.o
	$(CLINKER) $(OPTFLAGS) -o group_example group_example.o $(FLIBS)

comm_split_example: comm_split_example.o
	$(CLINKER) $(OPTFLAGS) -o comm_split_example comm_split_example.o $(FLIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(F77) $(FFLAGS) -c $*.f
.f90.o:
	$(F90) $(FFLAGS) -c $*.f90
.SUFFIXES: .f90
