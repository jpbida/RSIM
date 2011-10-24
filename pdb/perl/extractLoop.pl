#!/usr/bin/perl
use strict;
### This program uses  the pdb_db directory                         ###
### and identifies all fragments that are contained in a helix      ###

### Criteria continuous residues numbers on each side ###
#
# Final Output Format  #

## pdb_id,res1a,res2a,res3a,chainA,res1b,res2b,res3b,chainB,face1a...,face1b..

my $pdb_dir=$ARGV[0];
my @files=<$pdb_dir/*_rna.ent_bps.txt>;

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
my $map_file=$file;
$map_file=~s/rna.ent/map.txt/g;
$map_file=~s/_bps.txt//g;
my $out_file=$map_file;
$out_file=~s/map.txt/loop.txt/g;
open(OT,">$out_file");
### Make sure map file only contains unique entries ####
system("sort $file | uniq > tmp$start.txt");
system("mv tmp$start.txt $file");

my $tri_file=$file;
$tri_file=~s/rna.ent/tri.txt/g;
$tri_file=~s/_bps.txt//g;
my $tri_file2=$tri_file;
$tri_file2=~s/tri.txt/ftri.txt/g;
print "$pdb_id: $map_file, $tri_file, $file\n";

open(MP, "<$map_file");
open(BP, "<$file" );

### Build Map Hash ###
my @maps=<MP>;
my @hashmap;
foreach my $map (@maps){
my $hash;
chomp($map);
$map=~s/ //g;
my @array=split(/,/,$map);
$hash->{"model"}=$array[0];
$hash->{"chain2"}=$array[1];
$hash->{"res2"}=$array[2];
$hash->{"res1"}=$array[4];
$hash->{"nuc"}=$array[5];
push(@hashmap,$hash);
}
close(MP);
### Map Bp file and build BP array ###
my @bps=<BP>;
my @bp_map;
foreach my $bp (@bps){
my $hash;
chomp($bp);
$bp=~s/ //g;
my @array=split(/,/,$bp);
$hash->{"model"}=($array[1]);
$hash->{"res1"}=$array[2];
$hash->{"res2"}=$array[3];
my $face1="N";
if($array[4]==0){$face1="W";}
if($array[4]==1){$face1="H";}
if($array[4]==2){$face1="S";}
my $face2="N";
if($array[5]==0){$face2="W";}
if($array[5]==1){$face2="H";}
if($array[5]==2){$face2="S";}
$hash->{"face1"}=$face1;
$hash->{"face2"}=$face2;
$hash->{"hydro"}=$array[6];
push(@bp_map,$hash);
}
## Re-map BP file ##
my @new_bp;
foreach my $bp_hash (@bp_map){
my $new_res1;
my $new_res2;
my $new_chain1;
my $new_chain2;
my $new_nuc1;
my $new_nuc2;
my $good1=0;
my $good2=0;
foreach my $map2 (@hashmap){
if($map2->{"res1"}==$bp_hash->{"res1"} && $map2->{"model"}==$bp_hash->{"model"}){
$new_res1=$map2->{"res2"};
$new_chain1=$map2->{"chain2"};
$new_nuc1=$map2->{"nuc"};
}
}
foreach my $map2 (@hashmap){
if($map2->{"res1"}==$bp_hash->{"res2"} && $map2->{"model"}==$bp_hash->{"model"}){
$new_res2=$map2->{"res2"};
$new_chain2=$map2->{"chain2"};
$new_nuc2=$map2->{"nuc"};
}
}
my $new_hash;
$new_hash->{"model"}=$bp_hash->{"model"};
$new_hash->{"res1"}=$new_res1;
$new_hash->{"res2"}=$new_res2;
$new_hash->{"chain1"}=$new_chain1;
$new_hash->{"chain2"}=$new_chain2;
$new_hash->{"face1"}=$bp_hash->{"face1"};
$new_hash->{"face2"}=$bp_hash->{"face2"};
$new_hash->{"hydro"}=$bp_hash->{"hydro"};
$new_hash->{"nuc1"}=$new_nuc1;
$new_hash->{"nuc2"}=$new_nuc2;

push(@new_bp,$new_hash);
}

### Parse BP information to identify helical regions ####
my @loops;
foreach my $bp1 (@new_bp){
## Residue Output ###
my $loop_len;
my $res;
my $chain;
my $model;
my $pdb;
### Determine if bp is in a loop ###
my $res1a = $bp1->{"res1"};
my $res1b = $bp1->{"res2"};
my $chain1= "$bp1->{'chain1'}";
my $chain2= "$bp1->{'chain2'}";
$model=  $bp1->{"model"};
if($res1a < $res1b){
if($chain1 eq $chain2 && ($res1b-$res1a) == 5){
my $hash;
$loop_len=5;
$res=$res1a;
$chain=$chain1;
$model = $model;
$pdb=$pdb_id;
print OT "$pdb_id,$model,$chain,$res,$loop_len\n";
}

if($chain1 eq $chain2 && ($res1b-$res1a) == 6){
$loop_len=6;
$res=$res1a;
$chain=$chain1;
$model = $model;
$pdb=$pdb_id;
print OT "$pdb_id,$model,$chain,$res,$loop_len\n";
}


}
}
close(OT);
close(BP);
close(MP);
#foreach $hash (@bp_map){
#my $model= $hash->{"model"};
#my $res1=$hash->{"res1"};
#my $res2=$hash->{"res2"};
#my $face1=$hash->{"face1"};
#my $face2=$hash->{"face2"};
#my $hydro=$hash->{"hydro"};
#print "$model,$res1,$res2,$face1,$face2,$hydro\n";
#}

### Add BP information to tri file ##



### write fuzzy.conf file ###
### execute pdb2fuzzy     ###
}
