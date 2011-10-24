#!/usr/bin/perl
my $ext=$ARGV[0];
my $ofile=$ARGV[1];

if(-e $ext){
open(FH,"<$ext");
open(OUT,">$ofile");

@lines=<FH>;

foreach $line (@lines){
### Get the sequence and generate all possible trinucleotides ###
## pdb,chainId,pos,seq ##
@sl=split(/,/,$line);
$seq=$sl[2];
$l=length($seq);
$i=0;
while($i<($l-3)){
$s1=substr($seq,$i,3);
$i++;
print OUT "$sl[0],$sl[1],$i,$s1\n"; 
}
}
close(OUT);
close(FH);
system("sort tris.dat -t',' -k4 > stris.dat")
}else{
printf "can not find input file: $ext\n";
}
