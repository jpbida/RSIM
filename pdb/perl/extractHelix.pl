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
$out_file=~s/map.txt/helix.txt/g;
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
my @helix;
foreach my $bp1 (@new_bp){
## Residue Output ###
my $res1a;
my $res2a=-1;
my $res3a=-1;
my $res1b;
my $res2b=-1;
my $res3b=-1;
my $face1a;
my $face2a;
my $face3a;
my $face1b;
my $face2b;
my $face3b;
my $model;
my $chain1;
my $chain2;
my $hydro1;
my $hydro2;
my $hydro3;
my $nuc1a;
my $nuc1b;
my $nuc2a;
my $nuc2b;
my $nuc3a;
my $nuc3b;
### determine if this bp is at the start of a trinucleotide helix ###
$nuc1a = $bp1->{"nuc1"};
$nuc1b = $bp1->{"nuc2"};
$res1a = $bp1->{"res1"};
$res1b = $bp1->{"res2"};
$face1a = $bp1->{"face1"};
$face1b = $bp1->{"face2"};
$chain1= "$bp1->{'chain1'}";
$chain2= "$bp1->{'chain2'}";
$model=  $bp1->{"model"};
$hydro1=$bp1->{"hydro"};
if($res1a < $res1b){
### Position 2 evaluation ###
my $loop1=0;
my $found=0;
while($loop1 <= $#new_bp){
my $model_t=$new_bp[$loop1]->{"model"};
my $nuc1_t=$new_bp[$loop1]->{"nuc1"};
my $nuc2_t=$new_bp[$loop1]->{"nuc2"};
my $res1_t=$new_bp[$loop1]->{"res1"};
my $res2_t=$new_bp[$loop1]->{"res2"};
my $chain1_t="$new_bp[$loop1]->{'chain1'}";
my $chain2_t="$new_bp[$loop1]->{'chain2'}";
my $face1_t=$new_bp[$loop1]->{"face1"};
my $face2_t=$new_bp[$loop1]->{"face2"};
my $hydro1_t=$new_bp[$loop1]->{"hydro"};

if($chain1 eq $chain1_t && $chain2 eq $chain2_t && $res1a == ($res1_t-1) && $res1b==($res2_t+1) && $model==$model_t){
$loop1=$#new_bp+1;
$nuc2a =$nuc1_t; 
$nuc2b =$nuc2_t; 
$res2a=$res1_t;
$res2b=$res2_t;
$face2a=$face1_t;
$face2b=$face2_t;
$hydro2=$hydro1_t;
$found=1;
}

$loop1++;
}
if($found==1){
$found=0;
$loop1=0;
while($loop1 <= $#new_bp){
my $model_t=$new_bp[$loop1]->{"model"};
my $nuc1_t=$new_bp[$loop1]->{"nuc1"};
my $nuc2_t=$new_bp[$loop1]->{"nuc2"};
my $res1_t=$new_bp[$loop1]->{"res1"};
my $res2_t=$new_bp[$loop1]->{"res2"};
my $chain1_t="$new_bp[$loop1]->{'chain1'}";
my $chain2_t="$new_bp[$loop1]->{'chain2'}";
my $face1_t=$new_bp[$loop1]->{"face1"};
my $face2_t=$new_bp[$loop1]->{"face2"};
my $hydro2_t=$new_bp[$loop1]->{"hydro"};
if($chain1 eq $chain1_t && $chain2 eq $chain2_t && $res2a == ($res1_t-1) && $res2b==($res2_t+1) && $model==$model_t){
$loop1=$#new_bp+1;
$nuc3a=$nuc1_t;
$nuc3b=$nuc2_t;
$res3a=$res1_t;
$res3b=$res2_t;
$face3a=$face1_t;
$face3b=$face2_t;
$hydro3=$hydro2_t;
$found=1;
}

$loop1++;
}
if($found==1){
print OT "$pdb_id,$model,$chain1,$chain2,$res1a,$res1b,$face1a:$face1b,$res2a,$res2b,$face2a:$face2b,$res3a,$res3b,$face3a:$face3b,$face1a-$face1b:$face2a-$face2b:$face3a-$face3b,$hydro1,$hydro2,$hydro3,$nuc1a$nuc1b,$nuc2a$nuc2b,$nuc3a$nuc3b\n";


}


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
