package bida::pdb;

use 5.008006;
use strict;
use warnings;
use Carp;

use bida::pdb::model;
use bida::pdb::atom;
use bida::pdb::chain;
use bida::pdb::residue;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use bida::rna ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.


### Exported if requested ##
our @EXPORT_OK = qw(
		);
### Exported Automatically ##
our @EXPORT = qw(
	
);

our $VERSION = '0.01';

## Maintains a list of files loaded
my @filelist;

### Store a vector of all files loaded ##
sub new
{
	my $class=shift;
my $args={@_};
bless($args,$class);
my $ar=[];
	my $self = {
	models=>$ar,
	file=>$args->{file}
	};
	bless($self,$class);
## Create a hash array that 	
$self->_init;	
	return $self;
}

sub _init{
	my $self = shift;
	push @filelist,$self->file();
### Read in a pdb file to an atoms hash ##
my $fname=$self->file();
open(FH,"<$fname");
my @lines=<FH>;
close(FH);
### Find indices for each start and stop for models ###
### Indentify start and stops for each chain  ###
my $start=-1;
my $end=0;
my $index=0;
my $firstmod=0;
foreach my $line(@lines){
if($firstmod==1){last;}
if($start<0){
## looking for the start of a new model
        if($line=~m/^MODEL /){
       		$start=$index;         
	}else{
		if($line=~m/^ENDMDL/){print "PDB file missing MODEL tag\n"; $start=0;}
	}
}else{
	if($line=~m/^ENDMDL/){
	print "Reading Model\n";
	my $new_model=bida::pdb::model->new($start,$index,\@lines);
	$end=$index;
	$self->add_model($new_model);
	$firstmod=1;	
	$start=-1;
	}else{
		if($line=~m/^MODEL /){
			print "PDB file missing ENDMDL\n";
		my $new_model=bida::pdb::model->new($start,$index,\@lines);
		$end=$index;
		$self->add_model($new_model);
		$firstmod=1;	
		}
	}

}
$index++;
}
## If only one model existed and didn't have tags
if($end==0){
		my $new_model=bida::pdb::model->new($start,$index,\@lines);
		$self->add_model($new_model);
}

}

# Get/Set attributes in the class
sub _gen_accessor {
		my $aname = shift;
		eval "sub $aname { defined \$_[1] ?
			\$_[0]->{$aname} = \$_[1] :
			\$_[0]->{$aname} }";
}
_gen_accessor($_) for qw(file models);

sub add_model {
my $self=shift;
my $new_model=shift;
if(defined $new_model){ push @{$self->{models}},$new_model;}else{
print "No model was given to add to pdb\n";
}
}




 
sub res_stats{
my $self=shift;
### hash containing res_name=>A,counts,atoms(hash atom_name,counts) ###
my $cur_counts=shift;
my $model=@{$self->{models}}[0];
	foreach my $chain (@{$model->{chains}}){
		foreach my $residue (@{$chain->{residues}}){
				my $add=0;
				my $c=0;
				foreach my $cur (@{$cur_counts}){
					if($cur->{res_name} eq $residue->{name}){
						@{$cur_counts}[$c]->{count}=@{$cur_counts}[$c]->{count}+1;
						$add=1;
						foreach my $atom (@{$residue->{atoms}}){
							my $add1=0;
							my $c1=0;
							foreach my $curatom (@{$cur_counts->[$c]->{atoms}}){
								if($curatom->{name} eq $atom->{name}){
									$cur_counts->[$c]->{atoms}->[$c1]->{count}=$cur_counts->[$c]->{atoms}->[$c1]->{count} +1;	
								$add1=1;
								}
								$c1++;
							}
							if($add1==0){
							## Add atom to the cur_counts hash	
							my $ahash1=();
							my $anem=$atom->{name};	
							$ahash1->{name}=$atom->{name};
							print "Added $anem\n";
								$ahash1->{count}=1;
							push @{$cur_counts->[$c]->{atoms}}, $ahash1;
							}	
						}
					}			
				$c++;
				}
				if($add==0){
					my $hash;
					$hash->{res_name}=$residue->{name};
					$hash->{count}=1;
					my @ar;
					foreach my $natom(@{$residue->{atoms}}){
						my $ahash;
						$ahash->{name}=$natom->{name};
						$ahash->{count}=1;
						push @ar,$ahash;
					}
					$hash->{atoms}=\@ar;
				push @{$cur_counts},$hash;	
				}	
			}
		}
}

sub print{
my $self=shift;
foreach my $model(@{$self->{models}}){
print "MODEL  \n";
$model->print;
print "ENDMDL  \n";
}
}


1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

bida::pdb - Perl extension for blah blah blah

=head1 SYNOPSIS

  use bida::rna;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for bida::rna, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

John Bida, E<lt>bida@apple.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by John Bida

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.6 or,
at your option, any later version of Perl 5 you may have available.


=cut
