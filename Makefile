
# ==================================
# define our source and object files
# ==================================

SOURCES=Spanner.cpp SpanDet.cpp RunControlParameterFile.cpp Function-Generic.cpp Function-Sequence.cpp Histo.cpp MosaikAlignment.cpp PairedData.cpp headerSpan.cpp cluster.cpp steps.cpp DepthCnvDet.cpp BedFile.cpp  SHA1.cpp
OBJECTS=$(SOURCES:.cpp=.o)

CSOURCES=fastlz.c
COBJECTS=$(CSOURCES:.c=.o)

# ================
# compiler options
# ================
RE2_ROOT=<edit this to be the area above the RE2 install area on your system>/re2
#RE2_ROOT=/usr/local/re2
BAMTOOLS_ROOT=<edit this to be the area above the RE2 install area on your system>/bamtools
#BAMTOOLS_ROOT=/usr/local/bamtools

CPPFLAGS=-Wall -O2 -march=nocona -I$(RE2_ROOT) -I$(BAMTOOLS_ROOT)/include -std=gnu++0x
 
LDFLAGS=-Wl,-s -static

PROGRAM=Spanner
LIBS=-L$(RE2_ROOT)/obj -L$(BAMTOOLS_ROOT)/lib  -lre2 -lbamtools -lz -lpthread

all: $(PROGRAM)
 
$(PROGRAM): $(OBJECTS) $(COBJECTS)
	@echo "- linking" $(PROGRAM)
	@$(CXX) $(LDFLAGS) $(FLAGS) -o $@ $^ $(LIBS) 

.PHONY: clean

clean:
	rm -f *.o $(PROGRAM) *~
