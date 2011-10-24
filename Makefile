DIRS = src apps
.PHONY: all clean

all : ${DIRS}
	for DIR in ${DIRS}; do \
		make -C $${DIR} all; \
	done
clean : ${DIRS}
	-rm -f src/*.o
	-rm -f bin/* 
