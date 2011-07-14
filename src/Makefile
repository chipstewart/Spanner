# ==================================
# define our source and object files
# ==================================

SOURCES=Spanner.cpp SpanDet.cpp RunControlParameterFile.cpp Function-Sequence.cpp UtilityFunctions.cpp Histo.cpp MosaikAlignment.cpp PairedData.cpp headerSpan.cpp cluster.cpp steps.cpp DepthCnvDet.cpp BedFile.cpp  SHA1.cpp
OBJECTS=$(SOURCES:.cpp=.o)

#CSOURCES=fastlz.c bgzf.c
CSOURCES=fastlz.c
COBJECTS=$(CSOURCES:.c=.o)
 
# ================
# compiler options
# ================
RE2_ROOT =/home/radon00/stewart/cga2/tools/re2
BAMTOOLS_ROOT=/home/radon00/stewart/cga2/tools/bamtools
TCLAP_ROOT=/home/radon00/stewart/cga2/tools

#CPPFLAGS=-Wall -O2 -march=nocona -std=c++0x -I$(RE2DIR) -I$(BAMTOOLS_ROOT)/include
#CPPFLAGS=-Wall -O2 -march=nocona -I$(RE2_ROOT) -I$(BAMTOOLS_ROOT)/include -std=gnu++0x
CPPFLAGS=-Wall -O2 -march=nocona -I$(RE2_ROOT) -I$(BAMTOOLS_ROOT)/include -I$(TCLAP_ROOT) 
#CPPFLAGS=-Wall -g -march=nocona -I$(RE2_ROOT) -I$(BAMTOOLS_ROOT)/include -std=gnu++0x

#CPPFLAGS=-Wall -g
LDFLAGS=-Wl,-s -static
#LDFLAGS=-Wl
#LDFLAGS=-static
PROGRAM=Spanner
LIBS=-L/usr/lib64 -L$(RE2_ROOT)/obj -L$(BAMTOOLS_ROOT)/lib  -lre2 -lbamtools -lz -lpthread

all: $(PROGRAM)
 
$(PROGRAM): $(OBJECTS) $(COBJECTS)
	@echo "- linking" $(PROGRAM)
	@$(CXX) $(LDFLAGS) $(FLAGS) -o $@ $^ $(LIBS) 

.PHONY: clean

clean:
	rm -f *.o $(PROGRAM) *~
