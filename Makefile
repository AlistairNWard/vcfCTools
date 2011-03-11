#OBJ_DIR = ./
HEADERS = bed.h \
          split.h \
          stats.h \
          tool_annotate.h \
          tool_filter.h \
          tool_intersect.h \
          tool_merge.h \
          tool_stats.h \
          tool_union.h \
          tool_unique.h \
	  tool_validate.h \
          tools.h \
          vcf.h \
	  vcfCTools_tool.h
SOURCES = bed.cpp \
          split.cpp \
          stats.cpp \
          tool_annotate.cpp \
          tool_filter.cpp \
          tool_intersect.cpp \
          tool_merge.cpp \
          tool_stats.cpp \
          tool_union.cpp \
          tool_unique.cpp \
          tool_validate.cpp \
          tools.cpp \
          vcf.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BIN_SOURCES = vcfCTools.cpp

BINS = $(BIN_SOURCES:.cpp=)

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(BINS): $(BIN_SOURCES) $(OBJECTS)
	$(CXX) $(OBJECTS) $@.cpp -o $@ $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	rm -f $(BINS) $(OBJECTS) vcfCTools.o

.PHONY: clean all
