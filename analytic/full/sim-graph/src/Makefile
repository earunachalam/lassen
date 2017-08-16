include ~/makefile.$(OS)
	
CXXFLAGS=-Wall -O4 -ggdb -std=c++11 -march=native \
		 -I../../../../shared/graph/cpp \
		 $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS=-lm -larmadillo

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/sim-graph.x

all: $(EXEC)

$(EXEC): $(wildcard *.cpp $(wildcard *.tpp)) $(wildcard *.h)
	$(CXX) $(CXXFLAGS) $(wildcard *.cpp) $(LINKFLAGS) -o $(EXEC) $(LIBS) 

clean:
	rm -f *.o $(EXEC)

.phony: all clean
