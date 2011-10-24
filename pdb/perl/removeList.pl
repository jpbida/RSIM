#!/usr/bin/perl


### Removes results from particular PDBS from the pdb directory given ###

my $main=$ARGV[0];
my $new =$ARGV[1];
my $list=$ARGV[2];
my $frag=$ARGV[3];

open(FH,"<$list");
my @list=<FH>;
foreach $line (@list){
chomp($line);
$pdb_id=$line;
print "$pdb_id\n";
system("rm $main/$pdb_id".".ent.gz_ftri.txt");
system("rm $main/$pdb_id".".ent.gz_map.txt");
system("rm $main/$pdb_id".".ent.gz_rna.ent");
system("rm $main/$pdb_id".".ent.gz_rna.ent_bps.txt");
system("rm $main/$pdb_id".".ent.gz_tri.txt");
system("cp $main/$pdb_id".".ent.gz $new/");
system("rm $frag/$pdb_id"."*");
}
