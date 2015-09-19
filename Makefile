# Set ODEOPTION on command line to (-DRK_SUITE | -DARK_ODE | -DBOOST_ODEINT).

# CC=g++
# Compiler is to be set on command line.

ifeq ($(CC),g++)
	WALL=-Wall
endif

# TODO: These MPI, ARK ODE and Boost odeint paths are to be adjusted on BG/L or BG/P.
CFLAGS = $(ODEOPTION) -g $(WALL) -Iext/rksuite-1.0 -I/usr/include/mpich -I/usr/include/boost
$(info CFLAGS is set to "${CFLAGS}")

LIBRKSUITE = librksuite.a

LIBS = -lm -lmpich -lhdf5 -lhdf5_hl -lsundials_arkode -lsundials_nvecserial
EXE = coupledCellsModel

SRCS = $(shell ls src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

.PHONY: all
all: $(EXE)

$(EXE): $(OBJS) $(LIBRKSUITE)
	$(CC) $(CFLAGS) $(OBJS) lib/$(LIBRKSUITE) $(LIBS) -o $(EXE)

rksuite: $(LIBRKSUITE)

$(LIBRKSUITE): ext/rksuite-1.0/rksuite.o
	ar ru lib/$@ $^
	ranlib lib/$@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -fv $(OBJS)
	rm -fv ext/rksuite-1.0/*.o
	rm -fv lib/$(LIBRKSUITE)
	rm -fv $(EXE)
