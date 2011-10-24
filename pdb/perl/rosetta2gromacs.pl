#!/usr/bin/perl

#### This files converts rosetta output pdbs to gromacs pdbs #####
# Rapid identification of functional structures from selex selected aptamers ###
# We aim to replace traditional biochemical methods for determination of structural domains
# with bioinformatic techniques. This allows for classif

my $file=$ARGV[0];
open(FH,"<$file");
my @lines=<FH>;
foreach $line(@lines){
	chomp($line);
if($line=~/^ATOM/){
	my $atoms=getAtoms($line);
my $resname=$atoms->{'resName'};
my $name=$atoms->{'name'};
if($resName eq "  A" | $resName eq " rA" | $resName eq "  C" | $resName eq " rC" | $resName eq " rG" | $resName eq "  G" | $resName eq " rU" | $resName eq "  U" ){
### Change Atom names ###
if($resName=~m/A/){$atoms->{'resName'}=" RA";}
if($resName=~m/C/){$atoms->{'resName'}=" RC";}
if($resName=~m/G/){$atoms->{'resName'}=" RG";}
if($resName=~m/U/){$atoms->{'resName'}=" RU";}
if($name eq " OP1"){$atoms->{'name'}=" O1P";}
if($name eq " OP2"){$atoms->{'name'}=" O2P";}
if($name eq " H5'"){$atoms->{'name'}="H5'1";}
if($name eq "H5''"){$atoms->{'name'}="H5'2";}
if($name eq "HO2'"){$atoms->{'name'}="HO'2";}
if($name eq " H2'"){$atoms->{'name'}="H2'1";}
if($resName eq "  U" | $resName eq " rU"){
if($name eq "  O2"){$atoms->{'name'}="   O";}
}
$out=tprint($atoms);
	print "$out\n";
}
}
}
### Read in file ###
# rA->RA
# OP1 -> O1P
# OP2 -> O2P
# H5' -> H5'1
# H5''-> H5'2
#ATOM     34  OP1    O1P    amber94_45   -0.77600     2
#ATOM     35  OP2    O2P    amber94_45   -0.77600     3
#ATOM     55  H5'    H5'1    amber94_19    0.06790     6  
#ATOM     56  H5''   H5'2    amber94_19    0.06790     7
#ATOM     59  H2'    H2'1    amber94_19    0.09720    30      
#ATOM     60  HO2'   HO'2    amber94_25    0.41860    32     

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



