#!/usr/bin/perl 

#### Program reduces a trinuc to 6 points ###
### c1-c3 centers of back bone ##
### n1-n3 centers of NTs ###
### Takes a file name as the argument returns the points ####

my $dir = $ARGV[0];
my $op = $ARGV[1];
open(AP,">$op");
close(AP);
open(AP,">>$op");
my @files=<$dir/*>;
$id=0;
foreach $file (@files){
open(FH,"<$file");
my @lines=<FH>;
### Coordinates for N1 #####
my $n1x=99999;
my $n1y=99999;
my $n1z=99999;

### Coordinates for N2 #####
my $n2x=99999;
my $n2y=99999;
my $n2z=99999;

### Coordinates for N3 #####
my $n3x=99999;
my $n3y=99999;
my $n3z=99999;

### Coordinates for C1 #####
my $c1x=99999;
my $c1y=99999;
my $c1z=99999;
### Cordinates for C2 #####
my $c2x=99999;
my $c2y=99999;
my $c2z=99999;
### Cordinates for C3 #####
my $c3x=99999;
my $c3y=99999;
my $c3z=99999;

#### get the coordinates of the backbone ####

$c1=getElement(\@lines,0," C1'");
$c1x=$c1->{"x"};
$c1y=$c1->{"y"};
$c1z=$c1->{"z"};
$c2=getElement(\@lines,1," C1'");
$c2x=$c2->{"x"};
$c2y=$c2->{"y"};
$c2z=$c2->{"z"};
$c3=getElement(\@lines,2," C1'");
$c3x=$c3->{"x"};
$c3y=$c3->{"y"};
$c3z=$c3->{"z"};


##### Get Coordinates of NTs #####
### Get the three positions and resNames ###

### Use get element ### 
my @array=get3NT(\@lines);
my $pt1=$array[0];
$n1x=$pt1->{"x"};
$n1y=$pt1->{"y"};
$n1z=$pt1->{"z"};
my $pt2=$array[1];
$n2x=$pt2->{"x"};
$n2y=$pt2->{"y"};
$n2z=$pt2->{"z"};
my $pt3=$array[2];
$n3x=$pt3->{"x"};
$n3y=$pt3->{"y"};
$n3z=$pt3->{"z"};
#### Get output ###
$out="";
$out.="$id,c1,$c1x,$c1y,$c1z\n";
$out.="$id,c2,$c2x,$c2y,$c2z\n";
$out.="$id,c3,$c3x,$c3y,$c3z\n";
$out.="$id,n1,$n1x,$n1y,$n1z\n";
$out.="$id,n2,$n2x,$n2y,$n2z\n";
$out.="$id,n3,$n3x,$n3y,$n3z\n";
print AP $out;
#A=C6
#G=C6
#U=C4
#C=C4
$id++;
}
close(AP);



sub get3NT{
#### Returns the center of the three NTs in the given trinuc ####
#### Reurns and array of hashes ####
#$c1=getElement(\@lines,1," C3'");
my $outh="NA";
my $l=shift(@_);
my @lines = @$l;
my @arry=();
foreach $line (@lines){
$hash=getAtoms($line);
$ename=$hash->{"name"};
$resname=$hash->{"resName"};
#### 4 Cases ######
#A=C6
#G=C6
#U=C4
#C=C4
if($resname=~m/A/){
if($ename eq " C6 "){push(@arry,$hash);}
}else{
if($resname=~m/G/){
if($ename eq " C6 "){push(@arry,$hash);}
}else{
if($resname=~m/C/){
if($ename eq " C4 "){push(@arry,$hash);}
}else{
if($resname=~m/U/){
if($ename eq " C4 "){push(@arry,$hash);}
}else{
}
}
}
}

}
return @arry;
}

sub getElement{
### If you give it an atom name, atom lines, and occurance the function returns the line for that atom ###
#$c1=getElement(\@lines,1," C3'");
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
