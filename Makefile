#OBJ_DIR = ./
HEADERS = bed.h \
          bedStructure.h \
          Fasta.h \
          genotype_info.h \
          info.h \
          intersect.h \
          output.h \
          samples.h \
          SmithWatermanGotoh.h \
          split.h \
          stats.h \
          tool_annotate.h \
          tool_filter.h \
          tool_intersect.h \
          tool_merge.h \
          tool_stats.h \
          tools.h \
          trim_alleles.h \
          variant.h \
          vcf.h \
          vcf_aux.h \
	  vcfCTools_tool.h
#          tool_distributions.h \
	  tool_validate.h
SOURCES = bed.cpp \
          bedStructure.cpp \
          Fasta.cpp \
          genotype_info.cpp \
          info.cpp \
          intersect.cpp \
          output.cpp \
          samples.cpp \
          SmithWatermanGotoh.cpp \
          split.cpp \
          stats.cpp \
          tool_annotate.cpp \
          tool_filter.cpp \
          tool_intersect.cpp \
          tool_merge.cpp \
          tool_stats.cpp \
          tools.cpp \
          trim_alleles.cpp \
          variant.cpp \
          vcf.cpp \
          vcf_aux.cpp
#          tool_distributions.cpp \
          tool_validate.cpp \

OBJECTS= $(SOURCES:.cpp=.o)

BIN_SOURCES = vcfCTools.cpp

BINS = $(BIN_SOURCES:.cpp=)

all: $(OBJECTS) $(BINS)

CXX = g++ -lm
#CXX = g++ -g -lm
CXXFLAGS = -O3

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(BINS): $(BIN_SOURCES) $(OBJECTS)
	$(CXX) $(OBJECTS) $@.cpp -o $@ $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	rm -f $(BINS) $(OBJECTS) vcfCTools.o

.PHONY: clean all
