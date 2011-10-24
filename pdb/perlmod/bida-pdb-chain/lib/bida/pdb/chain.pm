package bida::pdb::chain;

use 5.008002;
use strict;
use warnings;
use bida::pdb::residue;
use bida::pdb::atom;


require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use bida::pdb::chain ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';
sub new
{
my $class=shift;
my $starts=shift;
my $stops=shift;
my $lines=shift;
my $name=shift;
my $type=shift;
my $ar=[];
        my $self = {
        residues=>$ar,
	name=>$name,
	type=>$type
        };
        bless($self,$class);
## Create a hash array that     
$self->_init($starts,$stops,$lines);
        return $self;
}

sub _init{
        my $self = shift;
        my $starts = shift;
        my $stops = shift;
        my $p3 = shift;
my @lines=@{$p3};
my $start=-1;
my $cur_resseq=0;
my $cur_resname="";
my $index=$starts;
### Indentify start and stops for each chain  ###
while($index < $stops){ 
	my $line=$lines[$index];
        if($index>=$starts && $index <$stops){
	if($line=~m/^ATOM  /){
                my $atom=bida::pdb::atom->new($line);
        if($start<0){
		$cur_resseq=$atom->{resSeq};
		$cur_resname=$atom->{resName};
		$start=$index;
	}else{
        if($cur_resseq ne $atom->{resSeq}){
		my $new_res=bida::pdb::residue->new($start,$index,\@lines,$cur_resname,"unknown");
		$self->add_residue($new_res);	
		$start=$index;
		$cur_resseq=$atom->{resSeq};
		$cur_resname=$atom->{resName};
	}   
	}
	}
        }
$index++;
}
if($start >0){
		my $new_res=bida::pdb::residue->new($start,$index,\@lines,$cur_resname,"unknown");
		$self->add_residue($new_res);	
		$start=-1;

}

}

sub add_residue {
my $self=shift;
my $new_residue=shift;
if(defined $new_residue){ push @{$self->{residues}},$new_residue;}else{
print "No residue was given to add to chain\n";
}
}

sub _gen_accessor {
                my $aname = shift;
                eval "sub $aname { defined \$_[1] ?
                        \$_[0]->{$aname} = \$_[1] :
                        \$_[0]->{$aname} }";
}
_gen_accessor($_) for qw(residues);

sub print{
my $self=shift;
foreach my $residue(@{$self->residues}){
$residue->print;
}

}

sub print_no_h{
my $self=shift;
foreach my $residue(@{$self->residues}){
$residue->print_no_h();
}

}
# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

bida::pdb::chain - Perl extension for blah blah blah

=head1 SYNOPSIS

  use bida::pdb::chain;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for bida::pdb::chain, created by h2xs. It looks like the
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

John Bida;bida.john@mayo.edu;14489813, E<lt>m003206@mayo.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by John Bida;bida.john@mayo.edu;14489813

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
