SRCS = $(shell ls src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

CC = mpixlcxx
CFLAGS = -g -Wall -Icontrib
RKS_LIB = librksuite.a

LIBS = -lhdf5 -lhdf5_hl -lsundials_arkode -lsundials_nvecserial
EXE = coupledCellsModel

.PHONY: all
all: $(EXE)

$(EXE): $(OBJS) $(RKS_LIB)
	$(CC) $(CFLAGS) $(OBJS) $(RKS_LIB) $(LIBS) -o $(EXE)

rksuite: $(RKS_LIB)

$(RKS_LIB): contrib/rksuite.o
	ar ru $@ $^
	ranlib $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -fv $(OBJS)
	rm -fv contrib/*.o
	rm -fv $(RKS_LIB)
	rm -fv $(EXE)
