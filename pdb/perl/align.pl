#!/usr/bin/perl
use Math::Trig':radial';
my $dir = $ARGV[0];
print $dir;
@files=<$dir/*.ent>;
### Given a trinuc file it algins all atoms to a common origin ###
foreach $file (@files){
$all="";
### Get C3' of first NT ###
print "$file\n";
open(FH,"<$file");
$file=~s/$dir\///g;
print "output:$file\n";
open(OUT,">a$file");
my @lines=<FH>;
$tx=0;
$ty=0;
$tz=0;
$c1x=0;
$c1y=0;
$c1z=0;
$c1=getElement(\@lines,1," C3'");
$stop=0;
if($c1 eq "NA"){
$stop=1;
}else{
$p1x=$c1->{"x"};
$p1y=$c1->{"y"};
$p1z=$c1->{"z"};
}
$c1=getElement(\@lines,1," C1'");
if($c1 eq "NA"){
$stop=1;
}else{
$p2x=$c1->{"x"};
$p2y=$c1->{"y"};
$p2z=$c1->{"z"};
}
if($c1 eq "NA"){
$stop=1;
}else{
$c1=getElement(\@lines,1," C4'");
$p3x=$c1->{"x"};
$p3y=$c1->{"y"};
$p3z=$c1->{"z"};
}
### Translate C3' to the origin ####
if($stop==0){
foreach $line (@lines){
my $hash=getAtoms($line);
$hash=alignStructs($hash,$p1x,$p1y,$p1z,$p2x,$p2y,$p2z,$p3x,$p3y,$p3z);
my $out=tprint($hash);
$all.="$out\n";
}
}else{
print "BAD FILE\n";
$all="";
}
print OUT $all;
close(OUT);
}
#################################################################################
# Subroutines 
#################################################################################

sub getElement{
### If you give it an atom name, atom lines, and occurance the function returns the line for that atom ###
my $outh="NA";
my $l=shift(@_);
	my @lines = @$l;
	my $occurance=shift(@_);
	my $ele = shift(@_);
$lc=0;
	foreach $line (@lines){
$hash=getAtoms($line);
$ename=$hash->{"name"};
if($ename eq $ele){
#print "$ename:$ele\n";
if($lc==$occurance){
$outh=$hash;
$lc++;
}else{
$lc++;
}
}
}
return $outh;
}

sub alignStructs{
#### Function that aligns structure based on three points ###
### Translate point 1 to origin ##
### Rotates point 2 onto positive x-axis ##
### Rotates point 3 into positive xy plane ###
	my $hash=shift(@_);
	my $p1x=shift(@_);
	my $p1y=shift(@_);
	my $p1z=shift(@_);
	
	my $p2x=shift(@_);
	my $p2y=shift(@_);
	my $p2z=shift(@_);
	
	my $p3x=shift(@_);
	my $p3y=shift(@_);
	my $p3z=shift(@_);
$p2x=$p2x-$p1x;
$p2y=$p2y-$p1y;
$p2z=$p2z-$p1z;
($r,$angy,$angz)=cartesian_to_spherical($p2x, $p2y, $p2z);
$test=();
$test->{"x"}=$p2x;
$test->{"y"}=$p2y;
$test->{"z"}=$p2z;
$test=rotateZ($test,($angy));
$test=rotateY($test,(-1*$angz));
$test=rotateY($test,(1.570796));
$extraAng1=0;
$extraAng2=0;
#print "C1\n";
#print $test->{"x"};
#print ",";
#print $test->{"y"};
#print ",";
#print $test->{"z"};
#print "\n";
$xval=$test->{"x"};
if($xval<0){
	$extraAng1=3.141593;
}
$test=();
$test->{"x"}=$p3x;
$test->{"y"}=$p3y;
$test->{"z"}=$p3z;
$test=translate($test,$p1x,$p1y,$p1z);
### Rotate C1'to the positive x axis ##
$test=rotateZ($test,($angy));
$test=rotateY($test,(-1*$angz));
$test=rotateY($test,(1.570796));
$test=rotateZ($test,($extraAng1));
$p3x=$test->{"x"};
$p3y=$test->{"y"};
$p3z=$test->{"z"};
$angz4=atan(($p3z/$p3y));
$test=rotateX($test,(-1*$angz4));
$yval=$test->{"y"};
if($yval<0){
	$extraAng2=3.141593;
}
$test=rotateX($test,($extraAng2));

$hash=translate($hash,$p1x,$p1y,$p1z);
### Rotate C1'to the positive x axis ##
$hash=rotateZ($hash,($angy));
$hash=rotateY($hash,(-1*$angz));
$hash=rotateY($hash,(1.570796));
$hash=rotateZ($hash,($extraAng1));
### Rotate C4' into X/Y Plane ####
$hash=rotateX($hash,(-1*$angz4));
$hash=rotateX($hash,($extraAng2));
return $hash;
}


### Rotate C4' into positive z/x plane by rotating around x axis ###  
sub tprint{
	my $hash=shift(@_);
$hash->{"x"}= sprintf("%.3f", $hash->{"x"});
$hash->{"y"}= sprintf("%.3f", $hash->{"y"});
$hash->{"z"}= sprintf("%.3f", $hash->{"z"});
$hash->{"x"}= sprintf("%8.8s", $hash->{"x"});
$hash->{"y"}= sprintf("%8.8s", $hash->{"y"});
$hash->{"z"}= sprintf("%8.8s", $hash->{"z"});
$hash->{"serial"}=sprintf("%5.5s",$hash->{"serial"});
$hash->{"name"}=sprintf("%4.4s",$hash->{"name"});
$hash->{"altLoc"}=sprintf("%1.1s",$hash->{"altLoc"});
$hash->{"resName"}=sprintf("%3.3s",$hash->{"resName"});
$hash->{"chainID"}=sprintf("%1.1s",$hash->{"chainID"});
$hash->{"resSeq"}=sprintf("%4.4s",$hash->{"resSeq"});
$hash->{"iCode"}=sprintf("%1.1s",$hash->{"iCode"});
$hash->{"occupancy"}=sprintf("%6.2f",$hash->{"occupancy"});
$hash->{"tempFactor"}=sprintf("%6.2f",$hash->{"tempFactor"});
$hash->{"element"}=sprintf("%2.2s",$hash->{"element"});
$hash->{"charge"}=sprintf("%2.2s",$hash->{"charge"});
	
my $out="ATOM  ".$hash->{"serial"}." ".$hash->{"name"}."".$hash->{"altLoc"}."".$hash->{"resName"}." ".$hash->{"chainID"}."".$hash->{"resSeq"}."".$hash->{"iCode"}."   ".$hash->{"x"}."".$hash->{"y"}."".$hash->{"z"}."".$hash->{"occupancy"}."".$hash->{"tempFactor"}."          ".$hash->{"element"}."".$hash->{"charge"};
return $out;
}

sub getAtoms{
	my $line=shift(@_);
	my $hash=();
if($line=~m/^ATOM/){
#COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
#  7 - 11        Integer       serial       Atom  serial number.
#  13 - 16        Atom          name         Atom name.
#  17             Character     altLoc       Alternate location indicator.
# 18 - 20        Residue name  resName      Residue name.
#  22             Character     chainID      Chain identifier.
#  23 - 26        Integer       resSeq       Residue sequence number.
#  27             AChar         iCode        Code for insertion of residues.
#  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#  55 - 60        Real(6.2)     occupancy    Occupancy.
#  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#  77 - 78        LString(2)    element      Element symbol, right-justified.
#  79 - 80        LString(2)    charge       Charge  on the atom.
$serial=substr($line,6,5);
$name=substr($line,12,4);
$altLoc=substr($line,16,1);
$resName=substr($line,17,3);
$chainID=substr($line,21,1);
$resSeq=substr($line,22,4);
$iCode=substr($line,26,1);
$x=substr($line,30,8);
$y=substr($line,38,8);
$z=substr($line,46,8);
$occupancy=substr($line,54,6);
$tempFactor=substr($line,60,6);
$element=substr($line,76,2);
$charge=substr($line,78,2);
$hash->{"serial"}=$serial;
$hash->{"name"}=$name;
$hash->{"altLoc"}=$altLoc;
$hash->{"resName"}=$resName;
$hash->{"chainID"}=$chainID;
$hash->{"resSeq"}=$resSeq;
$hash->{"iCode"}=$iCode;
$hash->{"x"}=$x;
$hash->{"y"}=$y;
$hash->{"z"}=$z;
$hash->{"occupancy"}=$occupancy;
$hash->{"tempFactor"}=$tempFactor;
$hash->{"element"}=$element;
$hash->{"charge"}=$charge;
}
return $hash;
}

sub translate{
my $hash=shift(@_);
my $tx=shift(@_);
my $ty=shift(@_);
my $tz=shift(@_);
$hash->{"x"}=$hash->{"x"}-$tx;
$hash->{"y"}=$hash->{"y"}-$ty;
$hash->{"z"}=$hash->{"z"}-$tz;
return $hash;
}

sub rotateY{
my $hash=shift(@_);
my $ang=shift(@_);
$z1=$hash->{"z"};
$x1=$hash->{"x"};
#	about Y
$hash->{"x"}=$x1*cos($ang)+$z1*sin($ang);
$hash->{"z"}=$z1*cos($ang)-$x1*sin($ang);
return $hash;
}
sub rotateZ{
my $hash=shift(@_);
my $ang=shift(@_);
$y1=$hash->{"y"};
$x1=$hash->{"x"};
#	about Y
$hash->{"x"}=$x1*cos($ang)+$y1*sin($ang);
$hash->{"y"}= -1*$x1*sin($ang)+$y1*cos($ang);
return $hash;
}

sub rotateX{
	use Math::Trig;
my $hash=shift(@_);
my $ang=shift(@_);
$z1=$hash->{"z"};
$y1=$hash->{"y"};
#	about 

$hash->{"y"}=$y1*cos($ang)-$z1*sin($ang);
$hash->{"z"}=$y1*sin($ang)+$z1*cos($ang);
return $hash;
}
