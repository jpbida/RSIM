CC=g++
CPPFLAGS= -O3 -ffast-math -I/usr/local/include/boost
 
.PHONY: all
all : $(patsubst %.cpp,%.o,$(wildcard *.cpp))
%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<  
