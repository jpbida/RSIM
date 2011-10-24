#!/usr/bin/perl
#
#   TestExample.pl
#
use strict;
use example;

my @intV = (0..3);
my @dblV = (0.0, 1.5, 1.0, 1.5);

print "Integer Array Values: ";
map { print "$_ ,"; } @intV;

print "\nAverage = ", example::average(\@intV), "\n";

print "\nTry passing in an array of doubles: \n";
eval {"Average = ", example::average(\@dblV), "\n"; };
print "\tError encountered: $@\n" if ($@);
#
print "Double Array Values: ";
map { print "$_ ,"; } @dblV;

my $ref = example::half(\@dblV);
print "\nHalved Values: ";
map { print "$_ ,"; } @$ref;
print "\n\n";
#
my $str = "Hello World!\n";
print example::TestMod($str); 

print "\nConstant PGSd_EOS_AM = ", $example::PGSd_EOS_AM, "\n";
print "Constant PGSd_EOS_PM = ", $example::PGSd_EOS_PM, "\n\n";

