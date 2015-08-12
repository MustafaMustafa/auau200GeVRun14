OBJS = StFastSimUtilities.o test.o
EXE = test

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS     += -lEG

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++ -m64
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< 

clean:
	rm -f *.o
