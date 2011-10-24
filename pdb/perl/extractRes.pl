#!/usr/bin/perl 
#### THis program retrieves all residues and tells us the types of atoms that appear ####
my @last5=();
my @last5id=();
my @last5raw=();
#### Required Atoms in Residues ###
my @all_nts_AG=(
" P  ",
" OP1",
" OP2",
" O5'",
" C5'",
" C4'",
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",
" C1'",
" C4 ",
" N9 "
);

my @all_nts_UC=(
" P  ",
" OP1",
" OP2",
" O5'",
" C5'",
" C4'",
" O4'",
" C3'",
" O3'",
" C2'",
" O2'",
" C1'",
" C2 ",
" N1 ");

my $i=0; 
while($i < $#all_nts_AG){
$meet_cons[$i++]=1;
}

### Directory containing all the pdbs ###
$dir = $ARGV[0];
$outpath = $ARGV[1];
print $dir;
@files=<$dir/*.gz>;
### unzip files and extract info ##
foreach $file (@files){
$out="file,model,chain_id,res_id1,res_id2,res_id3,res_id4,res_id5,res1,res2,res3,res4,res5\n";
open(SL,">$file"."_tri.txt");
print SL $out;
close(SL);
open(SL,">>$file"."_tri.txt");
my $resSeq_num=0;
my $out_pdb="";
my $out_map="";

my $cur_res_atoms="";
$i=0;
while($i < $#all_nts_AG){
$meet_cons[$i++]=1;
}
$start=0;
@last5=();
@last5id=();
@last5raw=();
$out="";
print "$file\n";
system("cp $file temp.gz");
system("gunzip temp.gz");
open(FH,"<temp");
### Get lines starting with ATOM #####
@lines=<FH>;
my $model=0;
my $seqlist="";
my $cur_res="";
my $cur_res_id="";
my $cur_id="";
my $cur_chain="";
my $min_res="";
foreach $line (@lines){
if($line=~m/^MODEL/){
$out_pdb.=$line;
$start=0;
$resSeq_num=0;
$cur_res="";
$cur_res_id="";
$cur_id="";
$cur_chain="";
$min_res="";
}
if($line=~m/^ENDMDL/){
$sum=0;
$model++;
foreach $icr (@meet_cons){
$sum=$sum+$icr;
}
## Yes = write residue without spaces to seqlist

if($sum==0){
$out_pdb.=$cur_res_atoms;
$out_map.="$model,$cur_chain,$cur_id,J,$resSeq2,$cur_res\n";
$resSeq_num++;
push(@last5,"$cur_res");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms="";
}else{
## No = write - to seqlist
push(@last5,"-");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms="";
}
if($#last5==4){
if($last5[0] eq "-" || $last5[1] eq "-" || $last5[2] eq "-" || $last5[3] eq "-" || $last5[4] eq "-"){}else{ 
### Make sure the residue numbers are continuous ###
my $diff1=$last5id[1]-$last5id[0];
my $diff2=$last5id[2]-$last5id[1];
my $diff3=$last5id[3]-$last5id[2];
my $diff4=$last5id[4]-$last5id[3];
if($diff1==1 && $diff2==1 && $diff3==1 && $diff4==1){
$file=~m/\/([\d\D]{1,4})\.ent\.gz/;
$f1=$1;
$cur_chain=~s/ //g;
my $id1=$last5id[1];
$id1=~s/ //g;
my $id2=$last5id[3];
$id2=~s/ //g;
$file_raw="$outpath/$f1"."_$model"."_$cur_chain"."_$id1"."_$id2.ent";
print "$file_raw\n";
$out2="$last5raw[0]$last5raw[1]$last5raw[2]$last5raw[3]$last5raw[4]";
$out="$file,$model,$chainID,$last5[0],$last5[1],$last5[2],$last5[3],$last5[4],$last5id[0],$last5id[1],$last5id[2],$last5id[3],$last5id[4],$file_raw\n";
print SL $out;
open(RL,">$file_raw");
print RL $out2;
close(RL);
}


}}
$cur_chain="";
$cur_res="";
$cur_res_id="";
$cur_id="";
$seqlist="";
#open(RO,">>rout.dat");
#print RO $out;
$i=0;
while($i < $#meet_cons){
$meet_cons[$i++]=1;
}



$resSeq_num=0;
$cur_chain=~s/ //g;
$cur_chain="";
$cur_res="";
$cur_res_id="";
$cur_id="";
$seqlist="";
$min_res="";


$out_pdb.=$line;
}
if($line=~m/^ATOM/){
my $resSeqT=substr($line,22,4);
my $iCodeT=substr($line,26,1);
$resSeqT=~s/ //g;
if($iCodeT ne " "){print "I hate PDB files! iCode is present: $resSeqT$iCodeT\n";}else{
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
$resSeq2=sprintf("%4.4s",$resSeq_num);
$chain2=sprintf("%1.1s","J");
my $resSeq_num2=$resSeq_num+1;

$resSeq3=sprintf("%4.4s",$resSeq_num2);

$x= sprintf("%8.8s", $x);
$y= sprintf("%8.8s", $y);
$z= sprintf("%8.8s", $z);
$serial=sprintf("%5.5s",$serial);
$name=sprintf("%4.4s",$name);
$altLoc=sprintf("%1.1s",$altLoc);
$resName=sprintf("%3.3s",$resName);
$chainID=sprintf("%1.1s",$chainID);
$resSeq=sprintf("%4.4s",$resSeq);
$iCode=sprintf("%1.1s",$iCode);
$occupancy=sprintf("%6.2f",$occupancy);
$tempFactor=sprintf("%6.2f",$tempFactor);
$element=sprintf("%2.2s",$element);
$charge=sprintf("%2.2s",$charge);
$line2="ATOM  $serial $name$altLoc$resName $chain2$resSeq2$iCode   $x$y$z$occupancy$tempFactor          $element$charge\n";
$line3="ATOM  $serial $name$altLoc$resName $chain2$resSeq3$iCode   $x$y$z$occupancy$tempFactor          $element$charge\n";
chomp($resName);
chomp($resSeq);

### Check to see if a given residue has all the right atoms  ##
if($start==0){
$cur_chain=$chainID;
$min_res=$resSeq;
$start=1;
#@last5=($resName);
#@last5id=($resSeq);
$cur_id=$resSeq;
$cur_res=$resName;
$cur_res_id=$resSeq;
}else{
if($#last5==4){
#### Print out a row in the database ####
$out="$file,$model,$chainID,$last5[0],$last5[1],$last5[2],$last5[3],$last5[4],$last5id[0],$last5id[1],$last5id[2],$last5id[3],$last5id[4]\n";
if($last5[0] eq "-" || $last5[1] eq "-" || $last5[2] eq "-" || $last5[3] eq "-" || $last5[4] eq "-"){}else{ 
my $diff1=$last5id[1]-$last5id[0];
my $diff2=$last5id[2]-$last5id[1];
my $diff3=$last5id[3]-$last5id[2];
my $diff4=$last5id[4]-$last5id[3];
if($diff1==1 && $diff2==1 && $diff3==1 && $diff4==1){
$file=~m/\/([\d\D]{1,4})\.ent\.gz/;
$f1=$1;
$cur_chain=~s/ //g;
my $id1=$last5id[1];
$id1=~s/ //g;
my $id2=$last5id[3];
$id2=~s/ //g;
$file_raw="$outpath/$f1"."_$model"."_$cur_chain"."_$id1"."_$id2.ent";
print "$file_raw\n";
$out2="$last5raw[0]$last5raw[1]$last5raw[2]$last5raw[3]$last5raw[4]";
$out="$file,$model,$chainID,$last5[0],$last5[1],$last5[2],$last5[3],$last5[4],$last5id[0],$last5id[1],$last5id[2],$last5id[3],$last5id[4],$file_raw\n";
print SL $out;
open(RL,">$file_raw");
print RL $out2;
close(RL);
}

}
my @last5tmp=@last5;
my @last5idtmp=@last5id;
my @last5rawtmp=@last5raw;
@last5=();
@last5id=();
@last5raw=();
$last5[0]=$last5tmp[1];
$last5[1]=$last5tmp[2];
$last5[2]=$last5tmp[3];
$last5[3]=$last5tmp[4];
$last5id[0]=$last5idtmp[1];
$last5id[1]=$last5idtmp[2];
$last5id[2]=$last5idtmp[3];
$last5id[3]=$last5idtmp[4];

$last5raw[0]=$last5rawtmp[1];
$last5raw[1]=$last5rawtmp[2];
$last5raw[2]=$last5rawtmp[3];
$last5raw[3]=$last5rawtmp[4];
}
}
if($cur_chain eq $chainID ){
if($cur_id eq $resSeq){
$cur_res_atoms.=$line2;
if($cur_res=~m/[AGag]/){
my $i=0;

foreach $atom (@all_nts_AG){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}

}else{

if($cur_res=~m/[UCuc]/){
my $i=0;
foreach $atom (@all_nts_UC){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}
}

}


}else{

### Did the residue meet all the constraints? ##
$sum=0;
foreach $icr (@meet_cons){
$sum=$sum+$icr;
}
## Yes = write residue without spaces to seqlist

if($sum==0){
$out_pdb.=$cur_res_atoms;
$out_map.="$model,$cur_chain,$cur_id,J,$resSeq2,$cur_res\n";
### Update with line 3 with ###
$resSeq_num++;
push(@last5,"$cur_res");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms=$line3;
}else{
## No = write - to seqlist
push(@last5,"-");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms=$line2;
}
# Move to next resName ##
$cur_res=$resName;
$cur_id=$resSeq;
# Reset $meet_cons 
my $i=0; 
while($i < $#meet_cons){
$meet_cons[$i++]=1;
}
if($cur_res=~m/[AGag]/){
my $i=0;

foreach $atom (@all_nts_AG){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}

}else{

if($cur_res=~m/[UCuc]/){
my $i=0;
foreach $atom (@all_nts_UC){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}
}

}


}
}else{
### Did the residue meet all the constraints? ##
$sum=0;
foreach $icr (@meet_cons){
$sum=$sum+$icr;
}
## Yes = write residue without spaces to seqlist

if($sum==0){
$out_pdb.=$cur_res_atoms;
$out_map.="$model,$cur_chain,$cur_id,J,$resSeq2,$cur_res\n";
### Update with line 3 with ###
$resSeq_num++;
push(@last5,"$cur_res");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms=$line3;
}else{
## No = write - to seqlist
push(@last5,"-");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms=$line2;
}

my $i=0; 
while($i < $#meet_cons){
$meet_cons[$i++]=1;
}
$cur_chain="$chainID";
$min_res=$resSeq;
$cur_res="$resName";
$cur_id="$resSeq";
$seqlist="";

if($cur_id eq $resSeq){

if($cur_res=~m/[AGag]/){
my $i=0;

foreach $atom (@all_nts_AG){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}

}else{

if($cur_res=~m/[UCuc]/){
my $i=0;
foreach $atom (@all_nts_UC){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}
}

}


}else{
### Did the residue meet all the constraints? ##
$sum=0;
foreach $icr (@meet_cons){
$sum=$sum+$icr;
}
## Yes = write residue without spaces to seqlist
if($sum==0){
$out_pdb.=$cur_res_atoms;
$out_map.="$model,$cur_chain,$cur_id,J,$resSeq2,$cur_res\n";
$cur_res_atoms=$line3;
$resSeq_num++;
push(@last5,"$cur_res");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
}else{
## No = write - to seqlist
push(@last5,"-");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms=$line2;
}
# Move to next resName ##
$cur_res=$resName;
$cur_id=$resSeq;

# Reset $meet_cons 
my $i=0; 
while($i < $#meet_cons){
$meet_cons[$i++]=1;
}

if($cur_res=~m/[AGag]/){
my $i=0;

foreach $atom (@all_nts_AG){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}

}else{

if($cur_res=~m/[UCuc]/){
my $i=0;
foreach $atom (@all_nts_UC){
if($atom eq $name){$meet_cons[$i++]=0;}
$i++;
}
}

}


}



}



}
}
### Create an R dataset ###
}
$sum=0;
foreach $icr (@meet_cons){
$sum=$sum+$icr;
}
## Yes = write residue without spaces to seqlist

if($sum==0){
$out_pdb.=$cur_res_atoms;
$out_map.="$model,$cur_chain,$cur_id,J,$resSeq2,$cur_res\n";
$resSeq_num++;
push(@last5,"$cur_res");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
$cur_res_atoms="";
}else{
## No = write - to seqlist
$cur_res_atoms="";
push(@last5,"-");
push(@last5id,"$cur_id");
push(@last5raw,"$cur_res_atoms");
}
if($#last5==4){
$out="$file,$model,$chainID,$last5[0],$last5[1],$last5[2],$last5[3],$last5[4],$last5id[0],$last5id[1],$last5id[2],$last5id[3],$last5id[4]\n";
if($last5[0] eq "-" || $last5[1] eq "-" || $last5[2] eq "-" || $last5[3] eq "-" || $last5[4] eq "-"){}else{ 
### Make sure the residue numbers are continuous ###
my $diff1=$last5id[1]-$last5id[0];
my $diff2=$last5id[2]-$last5id[1];
my $diff3=$last5id[3]-$last5id[2];
my $diff4=$last5id[4]-$last5id[3];
if($diff1==1 && $diff2==1 && $diff3==1 && $diff4==1){
$file=~m/\/([\d\D]{1,4})\.ent\.gz/;
$f1=$1;
$cur_chain=~s/ //g;
my $id1=$last5id[1];
$id1=~s/ //g;
my $id2=$last5id[3];
$id2=~s/ //g;
$file_raw="$outpath/$f1"."_$model"."_$cur_chain"."_$id1"."_$id2.ent";
print "$file_raw\n";
$out2="$last5raw[0]$last5raw[1]$last5raw[2]$last5raw[3]$last5raw[4]";
$out="$file,$model,$chainID,$last5[0],$last5[1],$last5[2],$last5[3],$last5[4],$last5id[0],$last5id[1],$last5id[2],$last5id[3],$last5id[4],$file_raw\n";
print SL $out;
open(RL,">$file_raw");
print RL $out2;
close(RL);
}




}}
$cur_chain="";
$cur_res="";
$cur_id="";
$seqlist="";
#open(RO,">>rout.dat");
#print RO $out;
while($i < $#meet_cons){
$meet_cons[$i++]=1;
}
close(FH);
open(FH2,">$file"."_rna.ent");
print FH2 $out_pdb;
close(FH2);
open(FH3,">$file"."_map.txt");
print FH3 $out_map;
close(FH3);

$out_pdb="";
system("rm temp");
close(SL);
}












