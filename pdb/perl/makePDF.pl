#!/usr/bin/perl
#### Creates a probability distribution for each trinucleotide ####
$dir=$ARGV[0];
open(FO,"<stris.dat");
@lines=<FO>;
$st="XXX";
mkdir $st;
$j=0;
foreach $line (@lines){
@ary=split(/,/,$line);
$lt=$ary[3];
$lt=~s/ //g;
chomp($lt);
print "$ary[0],$ary[1],$ary[2],$ary[3]\n";
### Get coordinates from pdb file ### 
system("cp $dir/pdb$ary[0]\.ent\.gz temp.gz");
system("gunzip temp.gz");
open(PB,"<temp");
@pdb=<PB>;
$chainNum=0;
$pos=0;
foreach $row (@pdb){
if($row=~m/^ATOM/){
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
$old=$resSeq;
$old2=$chainID;
$serial=substr($row,6,5);
$name=substr($row,12,4);
$altLoc=substr($row,16,1);
$resName=substr($row,17,3);
$chainID=substr($row,21,1);
$resSeq=substr($row,22,4);
$iCode=substr($row,26,1);
$x=substr($row,30,8);
$y=substr($row,38,8);
$z=substr($row,46,8);
$occupancy=substr($row,54,6);
$tempFactor=substr($row,60,6);
$element=substr($row,76,2);
$charge=substr($row,78,2);
if($chainID eq $ary[1]){
if($inseq==1){
if($pos<3){
	$out.="$row";
}else{
### Convert to normalized coordinates ####
if($st eq $lt){}else{
##write pdb and start again ###
mkdir $lt;
$j=0;
$st=$lt;
}

	open(FO,">$st\/pdb_$j.ent");
	print FO "$out\n";
	close(FO);
	$j++;

$inseq=0;
$pos=0;
$out="";
}
}else{
if($chainNum==$ary[2]){
	$inseq=1;
	$out="";
$pos=0;
}
}
if($old2!=$chainID){
	$chainNum=0;
	$pos=0;
}
if($resSeq!=$old){
$pos++;
print "$resSeq,$name\n";
$chainNum++;
}
}
}

}
close(PB);
system("rm temp");
### Align to first NT ###
	
}
