#!/usr/bin/perl

### This program uses  the pdb_db directory                         ###
### to calculated the base-pairing context of each fragment         ###

my $pdb_dir=$ARGV[0];
my $bin_path=$ARGV[1];
my $start=$ARGV[2];
my $stop=$ARGV[3];

my @files=<$pdb_dir/*_rna.ent>;
my @nfiles;

my $i=0;
foreach my $file(@files){
if($i >= $start && $i < $stop){
push(@nfiles,$file);
}
$i++;
}

print "Total Files $#nfiles\n";
foreach my $file (@nfiles){
$file=~m/\/([\d\D]{1,4})\.ent/;
$pdb_id=$1;
### write fuzzy.conf file ###
### execute pdb2fuzzy     ###
$outfile="$file"."_bps.txt";
my $conf="output_prefix=$outfile\nnative=$file\npdb_id=$pdb_id";
open(TM,">temp$start.conf");
print TM $conf;
close(TM);
system("$bin_path/pdb2fuzzy -i temp$start.conf");
}
