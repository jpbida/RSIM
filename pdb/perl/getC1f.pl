#!/usr/bin/perl 
$id=0;
my $file=$ARGV[0];
my $num=$ARGV[1];
open(FH,"<$file");
my @lines=<FH>;
	getElement(\@lines,$num," C1'");


sub get3NT_n2{
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
			if($ename eq " C4 "){push(@arry,$hash);}
		}else{
			if($resname=~m/G/){
				if($ename eq " C4 "){push(@arry,$hash);}
			}else{
				if($resname=~m/C/){
					if($ename eq " C2 "){push(@arry,$hash);}
				}else{
					if($resname=~m/U/){
						if($ename eq " C2 "){push(@arry,$hash);}
					}else{
					}
				}
			}
		}

	}
	return @arry;
}


sub get3NT_n2{
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
if($ename eq " C4 "){push(@arry,$hash);}
}else{
if($resname=~m/G/){
if($ename eq " C4 "){push(@arry,$hash);}
}else{
if($resname=~m/C/){
if($ename eq " C2 "){push(@arry,$hash);}
}else{
if($resname=~m/U/){
if($ename eq " C2 "){push(@arry,$hash);}
}else{
}
}
}
}

}
return @arry;
}

sub get3NT_n1{
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
if($ename eq " N9 "){push(@arry,$hash);}
}else{
if($resname=~m/G/){
if($ename eq " N9 "){push(@arry,$hash);}
}else{
if($resname=~m/C/){
if($ename eq " N1 "){push(@arry,$hash);}
}else{
if($resname=~m/U/){
if($ename eq " N1 "){push(@arry,$hash);}
}else{
}
}
}
}

}
return @arry;
}
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
$lti=0;
while($lc <= $occurance){
$line=$lines[$lti++];
if($lti==$#lines){
$lc=$occurance+1;
}
$hash=getAtoms($line);
$ename=$hash->{"name"};
if($ename eq $ele){
$lc++;
	my $x=$hash->{'x'};
	my $y=$hash->{'y'};
	my $z=$hash->{'z'};
	my $res=$hash->{'resName'};
	my $pos=$hash->{'resSeq'};
	print "$pos,$res,$x,$y,$z\n";
}
}
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

sub dist{
use Math::Trig;
use Math::VectorReal;
use constant PI    => 4 * atan2(1, 1);
my	$p1=shift(@_);
my	$p2=shift(@_);
my @v1=($p1->{'x'},$p1->{'y'},$p1->{'z'});
my @v2=($p2->{'x'},$p2->{'y'},$p2->{'z'});
#print "pos1: ";
foreach $cord(@v1){
#print "$cord:   ";
}
#print "\n";

#print "pos2: ";
foreach $cord(@v2){
#print "$cord:   ";
}
#print "\n";

my $dist=sqrt(($v1[0]-$v2[0])**2+($v1[1]-$v2[1])**2+($v1[2]-$v2[2])**2);

return $dist;
}

sub bond{
use Math::Trig;
use Math::VectorReal;
use constant PI    => 4 * atan2(1, 1);
my	$p1=shift(@_);
my	$p2=shift(@_);
my	$p3=shift(@_);
my @v1=($p1->{'x'},$p1->{'y'},$p1->{'z'});
my @v2=($p2->{'x'},$p2->{'y'},$p2->{'z'});
my @v3=($p3->{'x'},$p3->{'y'},$p3->{'z'});
#### Create vectors ####
my $a=vector(($v1[0]-$v2[0]),($v1[1]-$v2[1]),($v1[2]-$v2[2]));
my $b=vector(($v3[0]-$v2[0]),($v3[1]-$v2[1]),($v3[2]-$v2[2]));
### Normal vectors to the two planes ###
my $l1 = $a->length;
my $l2 = $b->length;
my $bond = 0;
if($l2!=0 & $l1!=0){
$bond = acos(($a.$b)/($l1 * $l2));
}
return $bond;
}

sub tor{
#### Calculate Torsion Angles ####
use Math::Trig;
use Math::VectorReal;
use constant PI    => 4 * atan2(1, 1);
my	$p1=shift(@_);
my	$p2=shift(@_);
my	$p3=shift(@_);
my	$p4=shift(@_);
my @v1=($p1->{'x'},$p1->{'y'},$p1->{'z'});
my @v2=($p2->{'x'},$p2->{'y'},$p2->{'z'});
my @v3=($p3->{'x'},$p3->{'y'},$p3->{'z'});
my @v4=($p4->{'x'},$p4->{'y'},$p4->{'z'});
#### Create vectors ####
my $a=vector(($v2[0]-$v1[0]),($v2[1]-$v1[1]),($v2[2]-$v1[2]));
my $b=vector(($v3[0]-$v2[0]),($v3[1]-$v2[1]),($v3[2]-$v2[2]));
my $c=vector(($v4[0]-$v3[0]),($v4[1]-$v3[1]),($v4[2]-$v3[2]));

### Normal vectors to the two planes ###
my $n1=$a x $b;
my $n2=$b x $c;
my $l1 = $n1->length;
my $l2 = $n2->length;
my $tor = 0;
if($l2!=0 & $l1!=0){
$tor = acos(($n1.$n2)/($l1 * $l2));
#### Determining the sign of the angle ####
my $l3=$c->length;
my $sv = acos(($n1.$c)/($l1*$l3));
$m=PI;
$m=$m/2;
my $sign=0;
if($sv <= $m){
$sign = 1;
}else{
$sign = -1;
}
$tor = abs($tor)*$sign;
}
return $tor;

}



