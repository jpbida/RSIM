#!/usr/bin/perl

use bida::pdb;

### Read in a filelist and loop through all pdb files ###
my $pdb_id=$ARGV[0];
my $chain_id=$ARGV[1];
### Assume there aren't more than 200 residue types present in the pdb ###
	$pdb=bida::pdb->new(file=>"$pdb_id");
	foreach my $chain(@{$pdb->{models}->[0]->{chains}}){
		my $chain_name=$chain->{name};
	if($chain_name=~m/$chain_id/){	
		$chain->print_no_h();
	last;	
	}
}
