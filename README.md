# Makefile for compiling the N-Body simulation

CXX = g++
CXXFLAGS = -O2 -std=c++11
TARGET = nbody
SRC = n-body.cpp
OBJ = $(SRC:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ)

Run *make* to compile the simulation.

Executation:
Run *./nbody <num_partiles> <dt> <iterations> <output_interval>* to start the simulation.

Run *make clean* to remove complied files.
