#!/usr/bin/perl

#### Takes the information from bps2/ directory and renames the files in the AAA AAU ... directories ####
## ss_fname - single stranded RNA
## d123_fname - double stranded RNA at all positions
## d1_fname - double stranded RNA at position 1
## d12_fname - double stranded RNA at positions 1 2
## d13_fname - double stranded RNA at positions 1 3
## d23_fname - double stranded RNA at positions 2 3
## d2_fname - double stranded RNA at positions 2 
## d3_fname - double stranded RNA at positions 3 

### Read file that contains all directory names ####
open(FH,"<dirs.txt");
@dirs=<FH>;
close(FH);
foreach $dir (@dirs){
### get all files from the NT directory ###
	chomp($dir);
	@files=<$dir/*>;
foreach $file (@files){
### Read the file name ###
$p1=0;
$p2=0;
$p3=0;
$fname=$file;
$file=~s/$dir\///g;
$file=~s/\.ent//g;
@ary=split(/_/,$file);
### PDB, CHAIN, STARTPOS ###
print "$ary[0],$ary[1],$ary[2]\n";
$chainID=$ary[1];
$chainID=~s/ //g;
$startPos=$ary[2];
$startPos=~s/ //g;
$spos1=$startPos;
$spos2=$spos1+1;
$spos3=$spos1+2;
### Find the NTs for that trinuc in the bps2 directory 
if(-e "bps2\/$ary[0]\.bp"){
open(FO,"<bps2\/$ary[0]\.bp");
@lines=<FO>;
### Look for the start position, start+1, and start + 2 ####
foreach $line(@lines){
chomp($line);
@info=split(/,/,$line);
$chain2_2=$info[6];
$chain2_1=$info[5];
$chain2_1=~s/ //g;
$chain2_2=~s/ //g;
$pos2_1=$info[3];
$pos2_2=$info[4];
$pos2_1=~s/ //g;
$pos2_2=~s/ //g;
### must match chainID ####
if($chainID eq $chain2_1){
if($spos1 == $pos2_1){$p1=1;}
if($spos2 == $pos2_1){$p2=1;}
if($spos3 == $pos2_1){$p3=1;}
}else{
if($chainID eq $chain2_2){
if($spos1 == $pos2_2){$p1=1;}
if($spos2 == $pos2_2){$p2=1;}
if($spos3 == $pos2_2){$p3=1;}
}
}
}
close(FO);
}
## 000_fname - single stranded RNA
## 111_fname - double stranded RNA at all positions
## 100_fname - double stranded RNA at position 1
## 110_fname - double stranded RNA at positions 1 2
## 101_fname - double stranded RNA at positions 1 3
## 011_fname - double stranded RNA at positions 2 3
## 010_fname - double stranded RNA at positions 2 
## 001_fname - double stranded RNA at positions 3 
### rename the file ###
system("cp $fname pdbs\/$dir\_$p1$p2$p3\_$ary[0]\_$ary[1]_$ary[2]\.ent"); 
}
}

