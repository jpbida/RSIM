#!/usr/bin/perl 
### Directory containing all the pdbs ###
$dir = $ARGV[0];
@files=<$dir/*.gz>;
### unzip files and extract info ##
foreach $file (@files){
system("cp $file temp.gz");
system("gunzip temp.gz");
open(FH,"<temp");
### Get lines starting with ATOM #####
@lines=<FH>;
my $i=0;
while($i<$#lines){
my $line=$lines[$i];
$i++;
if($line=~m/^ATOM/){
my $resSeqT=substr($line,22,4);
my $iCodeT=substr($line,26,1);
$resSeqT=~s/ //g;
if($iCodeT ne " "){print "$file\n";$i=$#lines;}
}
}
close(FH);
system("rm temp");
}












