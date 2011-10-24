#!/usr/bin/perl
### Reads in a list of accession numbers and gets the zipped files from the pbd database

my $file=$ARGV[0];
my $ofile=$ARGV[1];
open(FH,"<$file");
my @acc_nums=<FH>;
foreach $acc_num (@acc_nums){ 
chomp($acc_num);
$acc_num=lc($acc_num);
$url="ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb$acc_num.ent.gz";
print "$url\n";
system("wget '$url' -O $ofile/$acc_num.ent.gz");
}
