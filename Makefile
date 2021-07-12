
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++11

main: main.cc *.hh
	$(CXX) $(CXXFLAGS) -o main main.cc

clean:
	rm -rf *.o main

