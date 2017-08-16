include ~/makefile.$(shell hostname)

CXXFLAGS= -Wall -O4 -g -ggdb -std=c++11 -march=native $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS=-lvoro++ -larmadillo -lfftw3 -lm

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/a.out

ALLOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
ALLOBJSDIR=$(addprefix $(OBJ)/, $(ALLOBJS))

all: $(EXEC)

$(EXEC): $(ALLOBJSDIR) $(wildcard *.h)
	$(CXX) $(ALLOBJSDIR) $(LINKFLAGS) -o $(EXEC) $(LIBS) 

$(OBJ)/main.o: main.cpp *.tpp *.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.gch $(EXEC)

print-%:
	echo $* = $($*)

.phony: all clean
