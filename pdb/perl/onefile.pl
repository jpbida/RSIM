#!/usr/bin/perl

###### Take the PDB's that had multiple models and just use the first model ####

### sall2.txt - contains pdb_ids with multiple models ####

open(FH,"<../sall2.txt");
@pdbs=<FH>;
foreach $pdb (@pdbs){
chomp($pdb);
system("cp \.\.\/split_nucs\/$pdb\_mod1\.ent \.\.\/nmod\/pdb$pdb\.ent");
system("gzip \.\.\/nmod\/pdb$pdb\.ent");
}
close(FH);
