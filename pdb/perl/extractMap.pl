#!/usr/bin/perl
use strict;
### This program uses  the pdb_db directory                         ###
### to calculated the base-pairing context of each fragment         ###

my $pdb_dir=$ARGV[0];
my @files=<$pdb_dir/*_rna.ent_bps.txt>;

print "Total PDB files: $#files\n";

foreach my $file (@files){
$file=~m/\/([\d\D]{1,4})\.ent/;
my $pdb_id=$1;
my $map_file=$file;
$map_file=~s/rna.ent/map.txt/g;
$map_file=~s/_bps.txt//g;

### Make sure map file only contains unique entries ####
system("sort $file | uniq > tmp.txt");
system("mv tmp.txt $file");


my $tri_file=$file;
$tri_file=~s/rna.ent/tri.txt/g;
$tri_file=~s/_bps.txt//g;
my $tri_file2=$tri_file;
$tri_file2=~s/tri.txt/ftri.txt/g;
print "$pdb_id: $map_file, $tri_file, $file\n";

open(MP, "<$map_file");
open(TR, "<$tri_file");
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
my $good1=0;
my $good2=0;
foreach my $map2 (@hashmap){
if($map2->{"res1"}==$bp_hash->{"res1"} && $map2->{"model"}==$bp_hash->{"model"}){
$new_res1=$map2->{"res2"};
$new_chain1=$map2->{"chain2"};
}
}
my $new_hash;
$new_hash->{"model"}=$bp_hash->{"model"};
$new_hash->{"res"}=$new_res1;
$new_hash->{"chain"}=$new_chain1;
$new_hash->{"face"}=$bp_hash->{"face1"}.".".$bp_hash->{"face2"};
push(@new_bp,$new_hash);
}

### Use the new_pb array to add base-pairing information to fragments ##
my @tris=<TR>;

#/pdb_db/1ldz.ent.gz,0,A,  C,  G,  A,  C,  C,   2,   3,   4,   5,   6,/1ldz_0_A_3_5.ent
## Skip header

my $start=0;
my @final_tri;
foreach my $tri (@tris){
if($start==0){
$start=1;
}else{
my $hash;
chomp($tri);
$tri=~s/ //g;
my @array=split(/,/,$tri);
$hash->{"file"}    =$array[0];
$hash->{"model"}   =$array[1];
$hash->{"chain_id"}=$array[2];
$hash->{"res1"}    =$array[3];
$hash->{"res2"}    =$array[4];
$hash->{"res3"}    =$array[5];
$hash->{"res4"}    =$array[6];
$hash->{"res5"}    =$array[7];
$hash->{"id1"}     =$array[8];
$hash->{"id2"}     =$array[9];
$hash->{"id3"}     =$array[10];
$hash->{"id4"}     =$array[11];
$hash->{"id5"}     =$array[12];
my $face1="-";
my $face2="-";
my $face3="-";
foreach my $map3 (@new_bp){
if($map3->{"model"}==$array[1] && $map3->{"chain"} eq $array[2] && $map3->{"res"}==$array[9]){$face1.=$map3->{"face"}.":";}
if($map3->{"model"}==$array[1] && $map3->{"chain"} eq $array[2] && $map3->{"res"}==$array[10]){$face2.=$map3->{"face"}.":";}
if($map3->{"model"}==$array[1] && $map3->{"chain"} eq $array[2] && $map3->{"res"}==$array[11]){$face3.=$map3->{"face"}.":";}
}
$hash->{"face1"}=$face1;
$hash->{"face2"}=$face2;
$hash->{"face3"}=$face3;
push(@final_tri,$hash);
}
}

close(TR);
open(TR,">$tri_file2");
foreach my $hash (@final_tri){
my $file=$hash->{"file"};
my $model=$hash->{"model"};
my $chain=$hash->{"chain_id"};
my $res1=$hash->{"res2"};
my $res2=$hash->{"res3"};
my $res3=$hash->{"res4"};
my $face1=$hash->{"face1"};
my $face2=$hash->{"face2"};
my $face3=$hash->{"face3"};
my $first_res=$hash->{"id2"};
$res1=~s/ //g;
$res2=~s/ //g;
$res3=~s/ //g;

print TR "$file,$model,$chain,$res1$res2$res3,$first_res,$face1,$face2,$face3\n";
 
}
close(TR);
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




}
