#!/usr/bin/perl
#### Creates a probability distribution for each trinucleotide ####
#### Directory of raw pdbs ####
$dir=$ARGV[0];

#### File containing pdb_id and 1st and 2nd strand NTs ####
## pdb_id,chain_id1,s1,s2,s3,chain_id2,d1,d2,d3 ##
open(FO,"<double.dat");
@lines=<FO>;
$j=0;
foreach $line (@lines){
chomp($line);
@ary=split(/,/,$line);
$chain_id1=$ary[1];
$s1=$ary[2];
$s2=$ary[3];
$s3=$ary[4];
$chain_id2=$ary[5];
$d1=$ary[6];
$d2=$ary[7];
$d3=$ary[8];
### Get coordinates from pdb file ### 
print "$dir/pdb$ary[0]\.ent\.gz\n";
system("cp $dir/pdb$ary[0]\.ent\.gz temp.gz");
system("gunzip temp.gz");
open(PB,"<temp");
@pdb=<PB>;
############################################################
$p1d=0;
$p2d=0;
$p3d=0;
my $out1s=getTriNuc($chain_id1,$s1,$s2,$s3,\@pdb);
my $tri=getTriNTs($chain_id1,$s1,$s2,$s3,\@pdb);
if($d1>=0){
$d1atom=getNuc($chain_id2,$d1,\@pdb);
$p1d=1;
}
if($d2>=0){
$d2atom=getNuc($chain_id2,$d2,\@pdb);
$p2d=1;
}
if($d3>=0){
$d3atom=getNuc($chain_id2,$d3,\@pdb);
$p3d=1;
}
### This contains both strands ####
### Print to the pdbs/AAA_111_pdb_chain_start ####

$fname="$tri\_$p1d$p2d$p3d\_$ary[0]_$chain_id1\_$s1\.ent";
open(OH,">dblpdbs/$fname");
print OH "$out1s$d1atom$d2atom$d3atom";
close(OH);
###################################
close(PB);
system("rm temp");
}


######################################################################
# Subroutines 
######################################################################

sub getTriNTs{
#### Function for returning all ATOMs in trinuc starting at a given resSeq in a given chain ###
my $chainID1=shift(@_);
my $pos1=shift(@_);
my $pos2=shift(@_);
my $pos3=shift(@_);
$pos1=sprintf("%4.4s",$pos1);
$pos2=sprintf("%4.4s",$pos2);
$pos3=sprintf("%4.4s",$pos3);
my $p1is=0;
my $p2is=0;
my $p3is=0;
my $l=shift(@_);
my @lines=@$l;
my $out="";
my $nt1="";
my $nt2="";
my $nt3="";
foreach my $row (@lines){
if($row=~m/^ATOM/){
my $serial=substr($row,6,5);
my $name=substr($row,12,4);
my $altLoc=substr($row,16,1);
my $resName=substr($row,17,3);
my $chainID=substr($row,21,1);
my $resSeq=substr($row,22,4);
my $iCode=substr($row,26,1);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
my $occupancy=substr($row,54,6);
my $tempFactor=substr($row,60,6);
my $element=substr($row,76,2);
my $charge=substr($row,78,2);
####### Extract ATOMS #####
#print "$chainID:$chainID1\n";
if($chainID eq $chainID1){
	if($resSeq eq $pos1){
$p1is=1;
$nt1=$resName;
}else{
if($resSeq eq $pos2){
$p2is=1;
$nt2=$resName;
}else{
if($resSeq eq $pos3){
$p3is=1;
$nt3=$resName;
}else{
if($p3is==1){
last;
}
}
}
}
}
	
}

}

$nt1=~s/ //g;
$nt2=~s/ //g;
$nt3=~s/ //g;
print "$nt1$nt2$nt3\n";
if($p1is==1 & $p2is==1 & $p3is==1){
$out="$nt1$nt2$nt3";
}else{
$out="NA";
}
return $out;
}



sub getTriNuc{
#### Function for returning all ATOMs in trinuc starting at a given resSeq in a given chain ###
my $chainID1=shift(@_);
#print "$chainID1\n";
my $pos1=shift(@_);
my $pos2=shift(@_);
my $pos3=shift(@_);
my $pos1=sprintf("%4.4s",$pos1);
my $pos2=sprintf("%4.4s",$pos2);
my $pos3=sprintf("%4.4s",$pos3);
my $p1is=0;
my $p2is=0;
my $p3is=0;
my $out="";
my $l=shift(@_);
my @lines=@$l;
foreach my $row (@lines){
if($row=~m/^ATOM/){
my $serial=substr($row,6,5);
my $name=substr($row,12,4);
my $altLoc=substr($row,16,1);
my $resName=substr($row,17,3);
my $chainID=substr($row,21,1);
my $resSeq=substr($row,22,4);
my $iCode=substr($row,26,1);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
my $occupancy=substr($row,54,6);
my $tempFactor=substr($row,60,6);
my $element=substr($row,76,2);
my $charge=substr($row,78,2);
####### Extract ATOMS #####
if($chainID eq $chainID1){
	if($resSeq eq $pos1){
#print "$resSeq:$pos1\n";
$p1is=1;
	$out.=$row;
}else{
if($resSeq eq $pos2){
$p2is=1;
$out.=$row;
}else{
if($resSeq eq $pos3){
$p3is=1;
$out.=$row;
}else{
if($p3is==1){
last;
}
}
}
}
}
	
}

}
if($p1is==1 & $p2is==1 & $p3is==1){
}else{
$out="NA";
}
return $out;
}


sub getNuc{
#### Function for returning all ATOMs in a given Nuc starting at a given resSeq in a given chain ###
my $chainID1=shift(@_);
#print "$chainID1\n";
my $pos1=shift(@_);
my $pos1=sprintf("%4.4s",$pos1);
my $p1is=0;
my $out="";
my $l=shift(@_);
my @lines=@$l;
foreach my $row (@lines){
if($row=~m/^ATOM/){
my $serial=substr($row,6,5);
my $name=substr($row,12,4);
my $altLoc=substr($row,16,1);
my $resName=substr($row,17,3);
my $chainID=substr($row,21,1);
my $resSeq=substr($row,22,4);
my $iCode=substr($row,26,1);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
my $occupancy=substr($row,54,6);
my $tempFactor=substr($row,60,6);
my $element=substr($row,76,2);
my $charge=substr($row,78,2);
####### Extract ATOMS #####
if($chainID eq $chainID1){
	if($resSeq eq $pos1){
#print "$resSeq:$pos1\n";
$p1is=1;
	$out.=$row;
}else{
if($p1is==1){last;}
}
}
}
}
if($p1is==1){}else{$out="NA";}
return $out;
}
