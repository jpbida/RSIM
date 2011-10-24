#!/usr/bin/perl

use bida::pdb;

### Read in a filelist and loop through all pdb files ###
$allpdbs=$ARGV[0];
print "$allpdbs\n";
open(FH,"<$allpdbs");
@filelist=<FH>;
### Assume there aren't more than 200 residue types present in the pdb ###
my @cur;
push @cur,"";
print "Length: $#cur\n";
foreach my $pdbfile (@filelist){
chomp($pdbfile);
	print "$pdbfile\n";
	$pdb=bida::pdb->new(file=>"$pdbfile");
	foreach my $chain(@{$pdb->{models}->[0]->{chains}}){
		my $chain_name=$chain->{name};
		my @chain_res;
		for($i=0; $i<200; $i++){push @chain_res, 0;}
		foreach my $residue(@{$chain->{residues}}){
			my $string = $residue->{name};
			$string.=$residue->OrderAtoms();
			my $eq=0;
			my $index=0;
			my $res_type=-1;
				while($index <= $#cur){
					if($cur[$index] eq $string){ $res_type=$index; $eq=1;$index=$#cur;}
					$index++;
				}
				if($eq==0 && $#cur < 200){push @cur,$string; print "Total $#cur\n";$res_type=$#cur;
open(RS,">>residues.txt");
print RS "$res_type $string\n";
close(RS);
}
$chain_res[$res_type]++;		
}
$scounts=join(" ",@chain_res);
$chain_out="$pdbfile $chain_name $scounts\n";
open(CS,">>chains.txt");
print CS "$chain_out";
close(CS);
}
	undef $pdb;
}
