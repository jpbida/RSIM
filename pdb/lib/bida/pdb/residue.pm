package bida::pdb::residue;

use 5.008002;
use strict;
use warnings;
use bida::pdb::atom;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use bida::pdb::residue ':all';
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
        atoms=>$ar,
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
my $index=$starts;
while($index<$stops){
my $line=$lines[$index];
        if($index>=$starts && $index <$stops){
        if($line=~m/^ATOM  /){
my $new_atom=bida::pdb::atom->new($line);
$self->add_atom($new_atom);
	}
}
$index++;
}
}

sub OrderAtoms{
my $self = shift;
my @anames;
	foreach my $key (@{$self->{atoms}}){
my $atom_name=$key->{name};
if($atom_name=~m/H/){}else{	
	push @anames,$key->{name};}	
	}
@anames=sort(@anames);
my $astring=join("",@anames);
return $astring;
}

sub get_type {
my $self=shift;
### Types - RNA, DNA, RNA-DNA, PROTEIN, OTHER ####



}


sub add_atom {
my $self=shift;
my $new_atom=shift;
if(defined $new_atom){ push @{$self->{atoms}},$new_atom;}else{
print "No atom was given to add to residue\n";
}
}

sub print{
my $self=shift;
foreach my $atom(@{$self->{atoms}}){
$atom->print;
}

}

sub print_no_h{
my $self=shift;
foreach my $atom(@{$self->{atoms}}){
$atom->print_no_h();
}

}

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

bida::pdb::residue - Perl extension for blah blah blah

=head1 SYNOPSIS

  use bida::pdb::residue;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for bida::pdb::residue, created by h2xs. It looks like the
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
