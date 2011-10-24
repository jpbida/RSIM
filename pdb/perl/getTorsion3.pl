#!/usr/bin/perl 
#### Program reduces a trinuc to 21 angles ###
#./getTorsion3.pl pdbs frag.db
my $dir = $ARGV[0];
$op="frag";
open(AP1,">$op\_db");
open(AP3,">$op\_nt");
open(AP2,">$op\_id");
open(AP4,">$op\_at");
close(AP1);
close(AP2);
close(AP3);
close(AP4);
open(AP1,">>$op\_db");
open(AP2,">>$op\_id");
open(AP3,">>$op\_nt");
open(AP4,">>$op\_at");
my @files=<$dir/*>;
$id=0;
foreach $file (@files){
print "$file\n";
open(FH,"<$file");
my @lines=<FH>;

#### Intialize arrays used ####
#T: len:ba:na1:na2:na3:nt

#1: len:  :   :   :na3:nt
#2: len:ba:na1:na2:na3:nt
#3:    :ba:na1:
### Get all backbone atoms ####
my @a1=getBBAtoms(" P  ",\@lines);
my @a2=getBBAtoms(" OP1",\@lines);
my @a3=getBBAtoms(" OP2",\@lines);
my @a4=getBBAtoms(" O5'",\@lines);
my @a5=getBBAtoms(" C5'",\@lines);
my @a6=getBBAtoms(" C4'",\@lines);
my @a7=getBBAtoms(" O4'",\@lines);
my @a8=getBBAtoms(" C3'",\@lines);
my @a9=getBBAtoms(" O3'",\@lines);
my @a10=getBBAtoms(" C2'",\@lines);
my @a11=getBBAtoms(" O2'",\@lines);
my @a12=getBBAtoms(" C1'",\@lines);
### Get 1-3 atoms ###
my $i=1;
my $out4="$id,$a1[$i]->{'x'},$a1[$i]->{'y'},$a1[$i]->{'z'},$a2[$i]->{'x'},$a2[$i]->{'y'},$a2[$i]->{'z'},$a3[$i]->{'x'},$a3[$i]->{'y'},$a3[$i]->{'z'},$a4[$i]->{'x'},$a4[$i]->{'y'},$a4[$i]->{'z'},$a5[$i]->{'x'},$a5[$i]->{'y'},$a5[$i]->{'z'},$a6[$i]->{'x'},$a6[$i]->{'y'},$a6[$i]->{'z'},$a7[$i]->{'x'},$a7[$i]->{'y'},$a7[$i]->{'z'},$a8[$i]->{'x'},$a8[$i]->{'y'},$a8[$i]->{'z'},$a9[$i]->{'x'},$a9[$i]->{'y'},$a9[$i]->{'z'},$a10[$i]->{'x'},$a10[$i]->{'y'},$a10[$i]->{'z'},$a11[$i]->{'x'},$a11[$i]->{'y'},$a11[$i]->{'z'},$a12[$i]->{'x'},$a12[$i]->{'y'},$a12[$i]->{'z'},";
$i=2;
$out4.="$a1[$i]->{'x'},$a1[$i]->{'y'},$a1[$i]->{'z'},$a2[$i]->{'x'},$a2[$i]->{'y'},$a2[$i]->{'z'},$a3[$i]->{'x'},$a3[$i]->{'y'},$a3[$i]->{'z'},$a4[$i]->{'x'},$a4[$i]->{'y'},$a4[$i]->{'z'},$a5[$i]->{'x'},$a5[$i]->{'y'},$a5[$i]->{'z'},$a6[$i]->{'x'},$a6[$i]->{'y'},$a6[$i]->{'z'},$a7[$i]->{'x'},$a7[$i]->{'y'},$a7[$i]->{'z'},$a8[$i]->{'x'},$a8[$i]->{'y'},$a8[$i]->{'z'},$a9[$i]->{'x'},$a9[$i]->{'y'},$a9[$i]->{'z'},$a10[$i]->{'x'},$a10[$i]->{'y'},$a10[$i]->{'z'},$a11[$i]->{'x'},$a11[$i]->{'y'},$a11[$i]->{'z'},$a12[$i]->{'x'},$a12[$i]->{'y'},$a12[$i]->{'z'},";
$i=3;
$out4.="$a1[$i]->{'x'},$a1[$i]->{'y'},$a1[$i]->{'z'},$a2[$i]->{'x'},$a2[$i]->{'y'},$a2[$i]->{'z'},$a3[$i]->{'x'},$a3[$i]->{'y'},$a3[$i]->{'z'},$a4[$i]->{'x'},$a4[$i]->{'y'},$a4[$i]->{'z'},$a5[$i]->{'x'},$a5[$i]->{'y'},$a5[$i]->{'z'},$a6[$i]->{'x'},$a6[$i]->{'y'},$a6[$i]->{'z'},$a7[$i]->{'x'},$a7[$i]->{'y'},$a7[$i]->{'z'},$a8[$i]->{'x'},$a8[$i]->{'y'},$a8[$i]->{'z'},$a9[$i]->{'x'},$a9[$i]->{'y'},$a9[$i]->{'z'},$a10[$i]->{'x'},$a10[$i]->{'y'},$a10[$i]->{'z'},$a11[$i]->{'x'},$a11[$i]->{'y'},$a11[$i]->{'z'},$a12[$i]->{'x'},$a12[$i]->{'y'},$a12[$i]->{'z'}\n";

#### get the coordinates of the backbone ####
my $n0=getElement(\@lines,0," C1'");
my $n2=getElement(\@lines,1," C1'");
my $n4=getElement(\@lines,2," C1'");
my $n6=getElement(\@lines,3," C1'");
my $n8=getElement(\@lines,4," C1'");
my $nt1=$n2->{'resName'};
my $nt2=$n4->{'resName'};
my $nt3=$n6->{'resName'};
### Use get element ### 
my @array1=get3NT_n1(\@lines);
my $n10=$array1[0];
my $n11=$array1[1];
my $n12=$array1[2];
my $n13=$array1[3];
my $n14=$array1[4];
my @array2=get3NT_n2(\@lines);
my $n20=$array2[0];
my $n21=$array2[1];
my $n22=$array2[2];
my $n23=$array2[3];
my $n24=$array2[4];

#### Calculate all torsion angles #####
#### Information for position 1 #######
my $len_1 = dist($n2,$n4); 
my $len_3 = dist($n6,$n8); 
my $len_0 = dist($n0,$n2); 
my $tor_1 = tor($n0,$n2,$n4,$n6);

#### Information for position 2 #######
my $len_2 = dist($n4,$n6);
my $tor_2 =  tor($n2,$n4,$n6,$n8);
my $bond1 = bond($n2, $n4, $n6); 
my $bond0 = bond($n0, $n2, $n4); 
my $bond2 = bond($n4, $n6, $n8); 
#### Information for position 3 #######
#### Align structure such that n4 is (0,0,0) n6 is (x,0,0) and n2 is (x,y,0) ###
$out1="$len_0\t$len_1\t$len_2\t$len_3\t$tor_1\t$tor_2\t$bond0\t$bond1\t$bond2\t$id\n";
my @ary=split(/_/,$file);
$out2="$id\t$ary[0]\t$ary[1]\t$ary[2]\t$ary[3]\t$ary[4]\n";
$id++;
my $x=$n2->{'x'};
my $y=$n2->{'y'};
my $z=$n2->{'z'};
my $bx=$n11->{'x'};
my $by=$n11->{'y'};
my $bz=$n11->{'z'};
my $b2x=$n21->{'x'};
my $b2y=$n21->{'y'};
my $b2z=$n21->{'z'};
$out3="$id,$x,$y,$z,$bx,$by,$bz,$b2x,$b2y,$b2z,$nt1,";
$x=$n4->{'x'};
$y=$n4->{'y'};
$z=$n4->{'z'};
$bx=$n12->{'x'};
$by=$n12->{'y'};
$bz=$n12->{'z'};
$b2x=$n22->{'x'};
$b2y=$n22->{'y'};
$b2z=$n22->{'z'};
$out3.="$x,$y,$z,$bx,$by,$bz,$b2x,$b2y,$b2z,$nt1,";
$x=$n6->{'x'};
$y=$n6->{'y'};
$z=$n6->{'z'};
$bx=$n13->{'x'};
$by=$n13->{'y'};
$bz=$n13->{'z'};
$b2x=$n23->{'x'};
$b2y=$n23->{'y'};
$b2z=$n23->{'z'};
$out3.="$x,$y,$z,$bx,$by,$bz,$b2x,$b2y,$b2z,$nt1\n";
print AP1 $out1;
print AP2 $out2;
print AP3 $out3;
print AP4 $out4;
}
close(AP1);
close(AP2);
close(AP3);
close(AP4);



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
sub getBBAtoms{
### Backbone atoms
my $outh="NA";
my $atomname=shift(@_);
my $l=shift(@_);
my @lines = @$l;
my @arry=();
foreach $line (@lines){
$hash=getAtoms($line);
$ename=$hash->{"name"};
$resname=$hash->{"resName"};
#### 4 Cases ######
if($ename eq $atomname){push(@arry,$hash);}
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



