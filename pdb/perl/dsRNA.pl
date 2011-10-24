#!/usr/bin/perl
#### Dectects dsRNA segments and classifies them with name.ds123ent,name.ds12,...name.ssent ###
## Where ds123 denotes the positions in the trinuc that are double stranded ###
### To dectect double stranded RNA we will use distances between 
#Foreach PDB file calculate all pairwise distances for G(N2,N1,O6), C(O2,N3,N4), A(N1,N6), U(N3,O4). Only using those atoms. 
#Extract all distances < 4A.  Distance file should look like
#PDB,NT1,POS1,NT2,POS2,ATOM1, ATOM2, Distance
my $dir=$ARGV[0];
open(TRI,"<sall3.dat");
@files=<TRI>;
foreach $fline (@files){
print "$fline\n";
chomp($fline);
$file1="$dir/pdb$fline\.ent\.gz";
system("cp $file1 temp.gz");
system("gunzip temp.gz");
my $file="temp";
open(FH,"<$file");
print "Reading in $file1\n";
my @lines=<FH>;
my $out="";
my $f1=$fline;
print "Parsing File $f1\n";
$out=parsePDB(\@lines);
@ary=split(/\n/,$out);
print "Using $#ary atoms\n";
### Calculate all pairwise distances in out ###
$j=1;
$pairs="";
foreach my $point (@ary)
{
chomp($point);
@ps=split(/,/,$point);
$left=($#ary)+1;
$m=$j;
$x1=$ps[3];
$y1=$ps[4];
$z1=$ps[5];
$chainID=$ps[6];
$z1=~s/ //g;
$y1=~s/ //g;
$x1=~s/ //g;
$po=$ps[1];
$po=~s/ //g;
chomp($z1);
while($m < $left){
	$p209=$ary[$m];
	chomp($p209);
	@ps2=split(/,/,$p209);
$x2=$ps2[3];
$y2=$ps2[4];
$z2=$ps2[5];
$z2=~s/ //g;
$y2=~s/ //g;
$x2=~s/ //g;
chomp($z2);
$x0=($x1-$x2)**2;
$y0=($y1-$y2)**2;
$z0=($z1-$z2)**2;
$dist=$z0+$x0+$y0;
$chainID2=$ps2[6];
$po0=$ps2[1];
$po0=~s/ //g;
$po1=$po0+1;
$po2=$po0-1;
if($dist<16){
$t1="$po";
$t2="$po0";
$t3="$po1";
$t4="$po2";
if($chainID eq $chainID2){
if($t1 eq $t2){}else{
	if($t1 eq $t3){}else{
	if($t1 eq $t4){}else{
$pairs.="$f1,$chainID,$ps[0],$ps[1],$ps[2],$chainID2,$ps2[0],$ps2[1],$ps2[2],$dist\n";
	}
	}
	}
}else{
$pairs.="$f1,$chainID,$ps[0],$ps[1],$ps[2],$chainID2,$ps2[0],$ps2[1],$ps2[2],$dist\n";
}

}
$m++;
	}
$j++;
}

### Convert Distances to basepairings ####
### Look only for watson-crick pairings ###
### G-C, A-U, G-U ###
### Write to file.bp ###
close(FH);
print "bps\/$f1\.bp\n";
open(OH,">bps\/$f1\.bp");
print OH $pairs;
close(OH);
system("rm temp");
}







sub parsePDB{
#### Function for returning all ATOMs in trinuc starting at a given resSeq in a given chain ###
my $l=shift(@_);
my @lines=@$l;
my $out = "";
foreach my $row (@lines){
	my $tmp=getBonds($row);
if($tmp eq "NABIDA"){
}else{
	$out.="$tmp\n";
}
}
return $out;
}

sub getBonds{
	$row=shift(@_);
	my $out="NABIDA";
my $chainID=substr($row,21,1);
my $name=substr($row,12,4);
my $resName=substr($row,17,3);
my $resSeq=substr($row,22,4);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
if($resName eq "  A"){
##########################
### Get Specific Atoms ###
##########################
if($name eq " N1 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}else{
if($name eq " N6 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
}
##########################
##########################

}else{
if($resName eq '  U'){
##########################
### Get Specific Atoms ###
##########################
	
if($name eq " N3 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}else{
if($name eq " O4 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
}
##########################
##########################

}else{
if($resName eq '  G'){
##########################
### Get Specific Atoms ###
##########################
if($name eq " N2 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
else{
if($name eq " N1 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}else{
if($name eq " O6 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
}
}
##########################
##########################
}else{
if($resName eq '  C'){
##########################
### Get Specific Atoms ###
##########################
if($name eq " O2 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
else{
if($name eq " N3 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}else{
if($name eq " N4 "){
$out="$resName,$resSeq,$name,$x,$y,$z,$chainID";
}
}
}	
}else{
$out="NABIDA";
}
}
}
}
return $out;
}
