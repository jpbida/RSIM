OBJECTS := $(wildcard ../src/*.o)
CC=g++
CPPFLAGS= -O3 -ffast-math -I/usr/include/boost -I../src 

.PHONY: all
all : $(patsubst %.cpp,%,$(wildcard *.cpp))
% : %.cpp $(OBJECTS)
	$(CC) $(CPPFLAGS) -o ../bin/$(patsubst %.cpp,%,$@) $@.cpp $(OBJECTS) 
