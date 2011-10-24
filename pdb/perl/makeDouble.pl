#!/usr/bin/perl
##### Makes the double.dat file for extracting paired RNA sequences ####
## pdb_id,chain_id1,s1,s2,s3,chain_id2,d1,d2,d3 ##
## Get all files from ./pdbs directory ##
@files=<pdbs/*.ent>;
$out="";
foreach $file (@files){
###CUG_111_1s9s_A_357.ent ###
chomp($file);
$file=~s/\.ent//g;
@ary=split(/_/,$file);
@pairs=split(//,$ary[1]);
$p1=$pairs[0];
$p2=$pairs[1];
$p3=$pairs[2];
if($p1==0 & $p2==0 & $p3==0){
}else{
$s2=$ary[4]+1;
$s3=$ary[4]+2;
$chain1=$ary[3];
### Open the bps2 file and get the positions of the coordinates ###
open(BP,"<bps2/$ary[2].bp");
#2a2e,  G,  C,1,334,A,A
@lines=<BP>;
$chain=-1;
if($p1==1){
$o1=getPair($ary[4],$chain1,\@lines);
($chain2_1,$p2_1)=split(/,/,$o1);
$chain=$chain2_1;
}
else{
$chain2_1= -1;
$p2_1=-1;
}
if($p2==1){
$o2=getPair($s2,$chain1,\@lines);
($chain2_2,$p2_2)=split(/,/,$o2);
if($chain == -1){
$chain=$chain2_2;
}else{
if($chain eq $chain2_2){}
else
{
print "Chains not equal $file: $chain2_2, $chain2_1\n";
}
}
}
else{
$chain2_2= -1;
$p2_2=-1;
}
if($p3==1){
$o3=getPair($s3,$chain1,\@lines);
($chain2_3,$p2_3)=split(/,/,$o3);
if($chain == -1){$chain=$chain2_3;}else{
if($chain eq $chain2_3){}else{print "Chains not equal $file: $chain2_3, $chain\n";}
}
}
else{
$chain2_3= -1;
$p2_3=-1;
}
### Need to check that the chain_id's are all the same ###
print "$ary[2],$ary[3],$ary[4],$s2,$s3,$chain,$p2_1,$p2_2,$p2_3\n";
}
}

sub getPair{
my $pos1=shift(@_);
my $chain1=shift(@_);
my $l=shift(@_);
my @lines=@$l;
foreach my $line(@lines){
chomp($line);
my @ary=split(/,/,$line);
if($ary[3]==$pos1 & $ary[5] eq $chain1){
$out="$ary[6],$ary[4]";
}else{
if($ary[4]==$pos1 & $ary[6] eq $chain1){
$out="$ary[5],$ary[3]";
}
}
}

#out = chain_id2,pos2 #
return $out;
}

