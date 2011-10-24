#!/usr/bin/perl
#### Creates a probability distribution for each trinucleotide ####
$dir=$ARGV[0];
open(FO,"<ntris.dat");
@lines=<FO>;
$j=0;
foreach $line (@lines){
@ary=split(/,/,$line);
$lt=$ary[1];
chomp($lt);
### Get coordinates from pdb file ### 
print "$dir/pdb$ary[0]\.ent\.gz\n";
system("cp $dir/pdb$ary[0]\.ent\.gz temp.gz");
system("gunzip temp.gz");
open(PB,"<temp");
@pdb=<PB>;
#### Getting minimum and maximum resSeq numbers ####
$minRes=9999;
$maxRes=0;
foreach $row (@pdb){
if($row=~m/^ATOM/){
$chainID=substr($row,21,1);
$resSeq=substr($row,22,4);
if($chainID eq $lt){
if($resSeq > $maxRes){
$maxRes=$resSeq;
}
if($resSeq < $minRes){
$minRes=$resSeq;
}
}

}

}
print "Minimum: $minRes, Maximum: $maxRes\n";
####
$start = $minRes;
$start=~s/ //g;
while($start <	$maxRes){
my $out1i=getTriNuc($lt,$start,($start+1),($start+2),\@pdb);
my $tri=getTriNTs($lt,$start,($start+1),($start+2),\@pdb);
print "$tri\n";
if($out1i eq "NA"){
print "No such resSeqs: $ary[0], chain $lt, start: $start\n";
}else{
print "$tri\n";
if (-d "$tri") {
}
else {
mkdir $tri;
}
open(FH,">$tri/$ary[0]_$lt\_$start.ent");
print FH $out1i;
close(FH);
}
$start++;
$tri="NA";
}


close(PB);
system("rm temp");
### Align to first NT ###
}


######################################################################
# Subroutines 
######################################################################

sub getTriNTs{
#### Function for returning all ATOMs in trinuc starting at a given resSeq in a given chain ###
my $chainID1=shift(@_);
my $pos1=shift(@_);
my $pos2=shift(@_);
my $pos3=shift(@_);
$pos1=sprintf("%4.4s",$pos1);
$pos2=sprintf("%4.4s",$pos2);
$pos3=sprintf("%4.4s",$pos3);
my $p1is=0;
my $p2is=0;
my $p3is=0;
my $l=shift(@_);
my @lines=@$l;
my $out="";
my $nt1="";
my $nt2="";
my $nt3="";
foreach my $row (@lines){
if($row=~m/^ATOM/){
my $serial=substr($row,6,5);
my $name=substr($row,12,4);
my $altLoc=substr($row,16,1);
my $resName=substr($row,17,3);
my $chainID=substr($row,21,1);
my $resSeq=substr($row,22,4);
my $iCode=substr($row,26,1);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
my $occupancy=substr($row,54,6);
my $tempFactor=substr($row,60,6);
my $element=substr($row,76,2);
my $charge=substr($row,78,2);
####### Extract ATOMS #####
#print "$chainID:$chainID1\n";
if($chainID eq $chainID1){
	if($resSeq eq $pos1){
$p1is=1;
$nt1=$resName;
}else{
if($resSeq eq $pos2){
$p2is=1;
$nt2=$resName;
}else{
if($resSeq eq $pos3){
$p3is=1;
$nt3=$resName;
}else{
if($p3is==1){
last;
}
}
}
}
}
	
}

}

$nt1=~s/ //g;
$nt2=~s/ //g;
$nt3=~s/ //g;
print "$nt1$nt2$nt3\n";
if($p1is==1 & $p2is==1 & $p3is==1){
$out="$nt1$nt2$nt3";
}else{
$out="NA";
}
return $out;
}



sub getTriNuc{
#### Function for returning all ATOMs in trinuc starting at a given resSeq in a given chain ###
my $chainID1=shift(@_);
#print "$chainID1\n";
my $pos1=shift(@_);
my $pos2=shift(@_);
my $pos3=shift(@_);
my $pos0=$pos1-1;
my $pos4=$pos3+1;
my $pos1=sprintf("%4.4s",$pos1);
my $pos2=sprintf("%4.4s",$pos2);
my $pos3=sprintf("%4.4s",$pos3);
my $pos4=sprintf("%4.4s",$pos4);
my $pos0=sprintf("%4.4s",$pos0);
my $p0is=0;
my $p4is=0;
my $p1is=0;
my $p2is=0;
my $p3is=0;
my $out="";
my $l=shift(@_);
my @lines=@$l;
foreach my $row (@lines){
if($row=~m/^ATOM/){
my $serial=substr($row,6,5);
my $name=substr($row,12,4);
my $altLoc=substr($row,16,1);
my $resName=substr($row,17,3);
my $chainID=substr($row,21,1);
my $resSeq=substr($row,22,4);
my $iCode=substr($row,26,1);
my $x=substr($row,30,8);
my $y=substr($row,38,8);
my $z=substr($row,46,8);
my $occupancy=substr($row,54,6);
my $tempFactor=substr($row,60,6);
my $element=substr($row,76,2);
my $charge=substr($row,78,2);
####### Extract ATOMS #####
if($chainID eq $chainID1){
	if($resSeq eq $pos0){$p0is=1;$out.=$row;}else{
		
	if($resSeq eq $pos4){$p4is=1;$out.=$row;}else{

	if($resSeq eq $pos1){$p1is=1;$out.=$row;}else{
if($resSeq eq $pos2){
$p2is=1;
$out.=$row;
}else{
if($resSeq eq $pos3){
$p3is=1;
$out.=$row;
}else{
if($p3is==1){
last;
}
}
}
}

}
}

}
	
}

}
if($p0is==1 & $p1is==1 & $p2is==1 & $p3is==1 & $p4is==1){
}else{
$out="NA";
}
return $out;
}


