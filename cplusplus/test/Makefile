## -*- Makefile -*-
##
## User: andrea
## Time: 4-gen-2017 16.49.53
## Makefile created by Oracle Developer Studio.
##
## This file is generated automatically.
##


#### Compiler and tool definitions shared by all build targets #####


CCC = g++
CXX = g++
CXXFLAGS = -std=c++11 -Wall -g


all : quadmesh_test
	

quadmesh_test : quadmesh_test.o quadmesh.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

quadmesh.o : ../src/quadmesh.cpp
	$(CXX) -c -o quadmesh.o ../src/quadmesh.cpp $(CXXFLAGS) 
	
quadmesh_test.o : quadmesh_test.cpp
	$(CXX) -c -o quadmesh_test.o quadmesh_test.cpp $(CXXFLAGS)

#%.o : %.cpp
#	$(CXX) -c -o %@ %< $(CXXFLAGS)

#### Clean target deletes all generated files ####
clean:
	rm -f quadmesh_test.o quadmesh.o quadmesh

