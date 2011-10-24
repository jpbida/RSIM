#!/bin/perl

my $dir = $ARGV[0];
my $prefix=$ARGV[1];
 
print "Extracting sequences from the given PDB files:"

./extractSeq.pl $dir $prefix
 
