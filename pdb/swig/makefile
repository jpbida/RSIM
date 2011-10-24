#
#       SWIG example make file for Linux 
#
all:
	swig -c++ -perl example.i
	g++ -c example.cpp example_wrap.cxx -fPIC\
        	`/usr/bin/perl -MExtUtils::Embed -e ccopts`
	g++ -shared example.o example_wrap.o -L /usr/lib/ \
                `/usr/bin/perl -MExtUtils::Embed -e ldopts`\
       		  -o example.so
        
clean:
	rm *.o *.cxx
        
purge:
	rm *.o *.so *.cxx *.pm

