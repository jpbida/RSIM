#!/usr/bin/perl
##### This program splits a given pdb file with mutliple models into several files containing only 1 model ####
$dir=$ARGV[0];
$ndir=$ARGV[1];
print "Directory: $dir\n";
@files=<$dir/*.gz>;
### unzip files and extract info ##
foreach $file (@files){
$file=~m/$dir\/pdb([\d\D]*)\.ent\.gz/;
$fname=$1;
system("cp $file temp.gz");
system("gunzip temp.gz");
open(FH,"<temp");
### Get lines starting with ATOM #####
@lines=<FH>;
$end=1;
$out="";
my $model=0;
foreach $line(@lines){
if($line=~m/^MODEL/){
$end=0;
$out=$line;
$line=~s/ //g;
$line=lc($line);
$line=~s/el//g;
chomp($line);
$nfname="$fname\_$line$model\.ent\n";
$model++;
print "$nfname\n";
}else{
if($line=~m/^ENDMDL/){
$end=1;
$out.=$line;
#open(OH,">../split_nucs/$nfname");
#print OH "$out";
#close(OH);
$out="";
}else{
if($end==0){
if($line=~m/^ATOM/){
$out.=$line;}
}
}
}
}
system("rm temp");
}
