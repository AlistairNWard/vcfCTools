#OBJ_DIR = ./
HEADERS = bed.h \
          Fasta.h \
          SmithWatermanGotoh.h \
          split.h \
          stats.h \
          tool_annotate.h \
          tool_filter.h \
          tool_intersect.h \
          tool_merge.h \
          tool_stats.h \
          tool_temp.h \
	  tool_validate.h \
          tools.h \
          vcf.h \
          vcf_aux.h \
	  vcfCTools_tool.h
SOURCES = bed.cpp \
          Fasta.cpp \
          SmithWatermanGotoh.cpp \
          split.cpp \
          stats.cpp \
          tool_annotate.cpp \
          tool_filter.cpp \
          tool_intersect.cpp \
          tool_merge.cpp \
          tool_stats.cpp \
          tool_temp.cpp \
          tool_validate.cpp \
          tools.cpp \
          vcf.cpp \
          vcf_aux.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BIN_SOURCES = vcfCTools.cpp

BINS = $(BIN_SOURCES:.cpp=)

all: $(OBJECTS) $(BINS)

CXX = g++ -lm
CXXFLAGS = -O3

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(BINS): $(BIN_SOURCES) $(OBJECTS)
	$(CXX) $(OBJECTS) $@.cpp -o $@ $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	rm -f $(BINS) $(OBJECTS) vcfCTools.o

.PHONY: clean all
