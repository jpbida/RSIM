#!/usr/bin/perl
use strict;
### This program uses  the pdb_db directory                         ###
### and identifies all fragments that are contained in a helix      ###

### Criteria continuous residues numbers on each side ###
#
# Final Output Format  #

## pdb_id,res1a,res2a,res3a,chainA,res1b,res2b,res3b,chainB,face1a...,face1b..

my $pdb_dir=$ARGV[0];
my @files=<$pdb_dir/*_tri.txt>;

print "Total PDB files: $#files\n";
my $start=$ARGV[1];
my $stop=$ARGV[2];
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
my $pdb_id=$1;
my $loop_file=$file;
my $helix_file=$file;
my $combo_file=$file;
$helix_file=~s/_tri/_helix/g;
$loop_file=~s/_tri/_loop/g;
$combo_file=~s/_tri/_combo/g;
### open helix.txt and loop.txt ###
### check if the frag is in a helix or a loop ##
open(HL,"<$helix_file");
open(LF,"<$loop_file");
open(FG,"<$file");
open(CM,">$combo_file");
my @loops=<LF>;
my @helix=<HL>;
my @frags=<FG>;  ## Skip first line
my @frags1=@frags[2..$#frags];
my @frags=@frags1;

foreach my $frag(@frags){
### Get Frag Information ###
chomp($frag);
my @array=split(/,/,$frag);
my $model=$array[1];
my $chain="$array[2]";
my $res1 =$array[4];
my $res2 =$array[5];
my $res3 =$array[6];
my $start=$array[9];
$res1=~s/ //g;
$res2=~s/ //g;
$res3=~s/ //g;
### Determine if it is in a helix ###
my $region_type=0; # 0 = other, 1=helix, 2=region
my $go=0;
while($go < $#helix){
my $hlix=$helix[$go];
$go=$go+1;
chomp($hlix);
my @helix_array=split(/,/,$hlix);
my $helix_model=$helix_array[1];
my $helix_chain="$helix_array[2]";
my $helix_start=$helix_array[4];
my $helix_type =$helix_array[13];
if($helix_start==$start && $chain eq $helix_chain && $model == $helix_model){
print CM "$pdb_id,$model,$chain,$start,$res1$res2$res3,1,$helix_type,0\n";
$go=$#helix;
$region_type=1;
}

}
if($region_type==0){
### Determine if it is in a loop ###
$go=0;
while($go < $#loops){
my $hlix=$loops[$go];
$go=$go+1;
chomp($hlix);
my @loop_array=split(/,/,$hlix);
my $loop_model=$loop_array[1];
my $loop_chain="$loop_array[2]";
my $loop_start=$loop_array[3];
my $loop_len =$loop_array[4];
if(($start-($loop_start-2)) < ($loop_len+2) && $chain eq $loop_chain && $model == $loop_model){
print CM "$pdb_id,$model,$chain,$start,$res1$res2$res3,2,---,$loop_len\n";
$go=$#loops;
$region_type=2;
}
}
if($region_type==0){
print CM "$pdb_id,$model,$chain,$start,$res1$res2$res3,0,---,0\n";
}
}

}
close(CM);
close(HL);
close(LF);
close(FG);
}
