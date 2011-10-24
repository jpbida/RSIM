#!/usr/bin/perl
use Data::Dumper;
$dir = $ARGV[0];
$prefix= $ARGV[1];
open(OH,">extract.dat");
print "Directory: $dir\n";
@files=<$dir/*.gz>;
### unzip files and extract info ##
foreach $file (@files){
$start=1;
$prints=0;
$seq="";
$out="";
$val=1;
$file=~m/$dir\/$prefix([\d\D]*)\.ent\.gz/;
$fname=$1;
system("cp $file temp.gz");
system("gunzip temp.gz");
open(FH,"<temp");
### Get lines starting with ATOM #####
@lines=<FH>;
foreach $line(@lines){
# COLUMNS      DATA TYPE       FIELD          DEFINITION
# -------------------------------------------------------------------
#  1 -  6      Record name     "SEQRES"
#  9 - 10      Integer         serNum         Serial number of the SEQRES record
#                                                  for the current chain. Starts at 1
#                                                  and increments by one each line.
#                                                  Reset to 1 for each chain.
#  12           Character       chainID        Chain identifier. This may be any
#                                                  single legal character, including a
#                                                  blank which is used if there is
#                                                  only one chain.
#  14 - 17      Integer         numRes         Number of residues in the chain.
#                                                  This value is repeated on every
#                                                  record.
#  20 - 22      Residue name    resName        Residue name.
#  24 - 26      Residue name    resName        Residue name.
#  28 - 30      Residue name    resName        Residue name.
#  32 - 34      Residue name    resName        Residue name.
#  36 - 38      Residue name    resName        Residue name.
#  40 - 42      Residue name    resName        Residue name.
#  44 - 46      Residue name    resName        Residue name.
#  48 - 50      Residue name    resName        Residue name.
#  52 - 54      Residue name    resName        Residue name.
#  56 - 58      Residue name    resName        Residue name.
#  60 - 62      Residue name    resName        Residue name.
#  64 - 66      Residue name    resName        Residue name.
#  68 - 70      Residue name    resName        Residue name.
#

$Record = substr($line,0,6);
if($Record=~m/SEQRES/){
@tary=seqRes($line);
#### Only Keep Sequences that are strictly RNA ####
if($start==1){
$start=0;}else{
if($tary[0]==1){
if($val==1){
print OH "$fname,$ary[1],$seq\n";
$seq="";
$prints=0;
}else{
### Reset the validation ###
$val=1;
}
$seq="";
}
}
@ary=seqRes($line);
$seq.="$ary[3]$ary[4]$ary[5]$ary[6]$ary[7]$ary[8]$ary[9]$ary[10]$ary[11]$ary[12]$ary[13]$ary[14]$ary[15]";
$prints=1;
$seq=~s/ //g;
$v=validSeq(\@ary);
if($v==0){
$val=0;
}
}else{
if($prints==1){
if($val==1){
print OH "$fname,$ary[1],$seq\n";
$seq="";
$prints=0;
}
}
}
}
system("rm temp");
}
close(OH);
sub seqRes{
$line=shift(@_);
@out=();
$serNum=substr($line,8,2);
$chainID=substr($line,11,1);
$numRes=substr($line,13,4);
$resName2=substr($line,19,3);
$resName3=substr($line,23,3);
$resName4=substr($line,27,3);
$resName5=substr($line,31,3);
$resName6=substr($line,35,3);
$resName7=substr($line,39,3);
$resName8=substr($line,43,3);
$resName9=substr($line,47,3);
$resName10=substr($line,51,3);
$resName11=substr($line,55,3);
$resName12=substr($line,59,3);
$resName13=substr($line,63,3);
$resName14=substr($line,67,3);
push(@out,$serNum);
push(@out,$chainID);
push(@out,$numRes);
push(@out,$resName2);
push(@out,$resName3);
push(@out,$resName4);
push(@out,$resName5);
push(@out,$resName6);
push(@out,$resName7);
push(@out,$resName8);
push(@out,$resName9);
push(@out,$resName10);
push(@out,$resName11);
push(@out,$resName12);
push(@out,$resName13);
push(@out,$resName14);
return @out;
}

sub validSeq{
$a=shift(@_);
@ar=@$a;
$valid=1;
$count=1;
foreach $ele (@ar){
	if($count > 3){
#print "$ele\n";
if($ele=~m/  A/){}else{
if($ele=~m/  C/){}else{
if($ele=~m/  G/){}else{
if($ele=~m/   /){}else{
if($ele=~m/  U/){}else{
	$valid=0;
}
}
}
}
}
}
$count++;
}
return $valid;
}
