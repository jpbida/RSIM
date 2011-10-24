#!/bin/bash

### Make all perl modules and install in the local library

cd /home/m003206/projects/RNA_Tools/pdb/perlmod/bida-pdb
perl Makefile.PL LIB=/home/m003206/projects/RNA_Tools/pdb/lib
make
make install

cd ../bida-pdb-model
perl Makefile.PL LIB=/home/m003206/projects/RNA_Tools/pdb/lib
make
make install

cd ../bida-pdb-chain
perl Makefile.PL LIB=/home/m003206/projects/RNA_Tools/pdb/lib
make
make install

cd ../bida-pdb-residue
perl Makefile.PL LIB=/home/m003206/projects/RNA_Tools/pdb/lib
make
make install

cd ../bida-pdb-atom
perl Makefile.PL LIB=/home/m003206/projects/RNA_Tools/pdb/lib
make
make install
