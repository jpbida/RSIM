#!/usr/bin/perl
use Math::Trig':radial';
use constant PI    => 4 * atan2(1, 1);
#### How are the structures stored? ########
# NT_id, 7 torsion angles, base_id
# Functions
# 	makePDB
# 		Takes torsion file and creates PDB
# 	checkSteric
# 		Checks all C1' positions for pairwise clashes
#
### Make structure an array of hashes 
	## ATOM =>
	## X=>
	## Y=>
	## Z=>
### Each NT is defined by N atoms ###
###
#### Reads in PDB and rotates torsion angles ####
# ./makeStruct.pl pdbs/
############################################################
#  Test of Torsion angle change                            #
############################################################
my $file = $ARGV[0];
#my $tor = $ARGV[1];
my $pos = $ARGV[2];
#### Load PDB into structure hash ####
open(FH,"<$file");
my @lines=<FH>;
#### Print out the inital hash ####
my @struct=();
foreach $line(@lines){
	my $points=getAtoms($line);
if($points->{'name'}=~m/O5'|C5'|C4'|C3'|O3'|P/){
#if($points->{'name'}=~m/P/){
	push(@struct,$points);
}
}


##### Can I use this to do the closed movements? ####
my $tor1=0;
my $max = 2*PI;
my $pos = 5;
$tor2=PI/10;
$im=0;
print "Loaded\n";
while($pos < 10){
	$tor1 = 0;
while($tor1 < $max){
print "MODEL $im\n";
$im++;
	changeTor(\@struct,$pos,$tor2);
	print "ENDMDL\n";
$tor1=$tor1+$tor2;
}
$pos = $pos+1;
}




sub changeTor{
### Array of hashes that define all the coordinates of a molecule ###
my $struct = shift(@_);
### Position of the 1st atom in the 4 that define the torsion angle #### 
my $pos1 = shift(@_);
my $tor = shift(@_);
my $pos2 = $pos1+1;
my $pos3 = $pos1+2;
my $pos4 = $pos1+3;

### p1 - p4 are hashes of coordinates ### 
my $p1 = $struct->[$pos1]; 
my $p2 = $struct->[$pos2]; 
my $p3 = $struct->[$pos3]; 
my $p4 = $struct->[$pos4];
my $p1x=$p1->{'x'};
my $p1y=$p1->{'y'};
my $p1z=$p1->{'z'};
my $p2x=$p2->{'x'};
my $p2y=$p2->{'y'};
my $p2z=$p2->{'z'};
my $p3x=$p3->{'x'};
my $p3y=$p3->{'y'};
my $p3z=$p3->{'z'};

### Align the entire molecule to the last three points ###
$i=0;
$all="";
foreach $point (@$struct){
$point=alignStructs($point,$p1x,$p1y,$p1z,$p2x,$p2y,$p2z,$p3x,$p3y,$p3z);
$x=$point->{'x'};
$y=$point->{'y'};
$z=$point->{'z'};
#print "Test2: $x,$y,$z\n";
### Rotate all points before pos1 around the x-axis ###
if($i<=$pos1){$point=rotateX($point,$tor);}
#my $out=tprint($point);
#$all.="$out\n";
$i++;
}
#print "$all\n";

}



sub makeStruct{
	my $hash=shift(@_);
### This program takes in at hash of torsion angles and produces an RNA structure ###
### Each NT has seven torsion angles ###
# alpha
# beta
# gamma
# delta
# epsilon
# zeta
# chi

### The input is an array of hashes ###
# The order of NT's in the array is the order they appear along 
# the chain
#

### Use standard bond lengths found in the five standards files 
# A.ent
# C.ent
# G.ent
# T.ent
# U.ent
# RNA.ent
#

#### For each angle #####

	

	
}

sub addBond{
### Given coordinates of A and B + torsion angle it outputs coordinates of C and D ###
my $hash = shift(@_);




}	

sub addSugar{
### Given C4', C3', Pucker, torsion, and NT it generates cordinates for O4', C1', C2' and base atoms ###



	
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
