CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra

SRCS = bam2iupac.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = bam2iupac

HTSLIB_DIR ?= /usr/local/lib

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) -lhts -Wl,-rpath,$(HTSLIB_DIR)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
