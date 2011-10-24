#!/usr/bin/perl
my $file=$ARGV[0];
my $prefix=$ARGV[1];
my $out_dir=$ARGV[2];
open(FH,"<$file");
my @acc_nums=<FH>;
foreach $acc_num (@acc_nums){ 
chomp($acc_num);
system("cp $prefix$acc_num.ent.gz $out_dir/$acc_num.ent.gz");
}
