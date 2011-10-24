#!/usr/bin/perl
use strict;
### This program uses  the pdb_db directory                         ###
### to calculated the base-pairing context of each fragment         ###

my $pdb_dir=$ARGV[0];
my @files=<$pdb_dir/*_rna.ent_bps.txt>;

print "Total PDB files: $#files\n";

### Run the r_script parseSS.r for each RNA molecule ###
foreach my $file (@files){
$file=~m/\/([\d\D]{1,4})\.ent/;
my $pdb_id=$1;
$cmd="R CMD BATCH '--args $pdb_dir $pdb_id'  parseSS.r"
print "Running $cmd\n";
system("$cmd");
}

