CC=g++
CPPFLAGS= -O3 -ffast-math -I /usr/local/include/boost
#CPPFLAGS= -O3 -ffast-math -I /home/m003206/lib/boost_1_42_0
#OPTIONS=g++ -O3 -ffast-math 
.PHONY: all
all : $(patsubst %.cpp,%.o,$(wildcard *.cpp))
%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<  
