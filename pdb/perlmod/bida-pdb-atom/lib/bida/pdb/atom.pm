package bida::pdb::atom;

use 5.008002;
use strict;
use warnings;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use bida::pdb::atom ':all';
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
### Takes in a line as an argument to create an atom ###
my $args=shift;
        my $self = {
	serial=>undef,
	name=>undef,
	altLoc=>undef,
	resName=>undef,
	chainID=>undef,
	resSeq=>undef,
	iCode=>undef,
	x=>undef,
	y=>undef,
	z=>undef,
	occupancy=>undef,
	tempFactor=>undef,
	element=>undef,
	charge=>undef
        };
        bless($self,$class);
## Create a hash array that     
$self->_init($args);
        return $self;
}
sub _init{
	my $self = shift;
	my $line = shift;
### Read in a pdb file to an atoms hash ##
$self->getAtoms($line);
}
sub _gen_accessor {
                my $aname = shift;
                eval "sub $aname { defined \$_[1] ?
                        \$_[0]->{$aname} = \$_[1] :
                        \$_[0]->{$aname} }";
}
_gen_accessor($_) for qw(serial name altLoc resName chainID resSeq iCode x y z occupancy tempFactor element charge);
sub print{
my $self=shift;
my $serial=$self->{"serial"};
my $name=$self->{"name"}			;
my $altLoc=$self->{"altLoc"};
my $resName=$self->{"resName"};
my $chainID=$self->{"chainID"};
my $resSeq=$self->{"resSeq"};
my $iCode=$self->{"iCode"};
my $x=$self->{"x"};
my $y=$self->{"y"};
my $z=$self->{"z"};
my $occupancy=$self->{"occupancy"};
my $tempFactor=$self->{"tempFactor"};
my $element=$self->{"element"};
my $charge=$self->{"charge"};
print "ATOM  $serial$name$altLoc$resName$chainID$resSeq$iCode$x$y$z$occupancy$tempFactor      $element$charge\n";
}

sub print_no_h{
my $self=shift;
my $serial=$self->{"serial"};
my $name=$self->{"name"}			;
my $altLoc=$self->{"altLoc"};
my $resName=$self->{"resName"};
my $chainID=$self->{"chainID"};
my $resSeq=$self->{"resSeq"};
my $iCode=$self->{"iCode"};
my $x=$self->{"x"};
my $y=$self->{"y"};
my $z=$self->{"z"};
my $occupancy=$self->{"occupancy"};
my $tempFactor=$self->{"tempFactor"};
my $element=$self->{"element"};
my $charge=$self->{"charge"};
if($name=~m/H/){}else{
print "ATOM  $serial$name$altLoc$resName$chainID$resSeq$iCode$x$y$z$occupancy$tempFactor      $element$charge\n";
}
}


sub getAtoms{
my $self=shift;
my $line=shift;
if($line=~m/^ATOM/){
#COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
#  7 - 11        Integer       serial       Atom  serial number.
#  13 - 16        Atom          name         Atom name.
#  17             Character     altLoc       Alternate location indicator.
# 18 - 20        Residue name  resName      Residue name.
#  22             Character     chainID      Chain identifier.
#  23 - 26        Integer       resSeq       Residue sequence number.
#  27             AChar         iCode        Code for insertion of residues.
#  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#  55 - 60        Real(6.2)     occupancy    Occupancy.
#  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#  77 - 78        LString(2)    element      Element symbol, right-justified.
#  79 - 80        LString(2)    charge       Charge  on the atom.
my $serial=substr($line,6,6);
my $name=substr($line,12,4);
my $altLoc=substr($line,16,1);
my $resName=substr($line,17,4);
my $chainID=substr($line,21,1);
my $resSeq=substr($line,22,4);
my $iCode=substr($line,26,4);
my $x=substr($line,30,8);
my $y=substr($line,38,8);
my $z=substr($line,46,8);
my $occupancy=substr($line,54,6);
my $tempFactor=substr($line,60,6);
my $element=substr($line,76,2);
my $charge=substr($line,78,2);
$self->{"serial"}=$serial;
$self->{"name"}=$name;
$self->{"altLoc"}=$altLoc;
$self->{"resName"}=$resName;
$self->{"chainID"}=$chainID;
$self->{"resSeq"}=$resSeq;
$self->{"iCode"}=$iCode;
$self->{"x"}=$x;
$self->{"y"}=$y;
$self->{"z"}=$z;
$self->{"occupancy"}=$occupancy;
$self->{"tempFactor"}=$tempFactor;
$self->{"element"}=$element;
$self->{"charge"}=$charge;
}
}


# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

bida::pdb::atom - Perl extension for blah blah blah

=head1 SYNOPSIS

  use bida::pdb::atom;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for bida::pdb::atom, created by h2xs. It looks like the
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
