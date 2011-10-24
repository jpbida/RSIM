# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
package example;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package examplec;
bootstrap example;
package example;
@EXPORT = qw( );

# ---------- BASE METHODS -------------

package example;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package example;

*average = *examplec::average;
*half = *examplec::half;
*TestMod = *examplec::TestMod;

############# Class : example::IntVector ##############

package example::IntVector;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( example );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = examplec::new_IntVector(@_);
    bless $self, $pkg if defined($self);
}

*size = *examplec::IntVector_size;
*empty = *examplec::IntVector_empty;
*clear = *examplec::IntVector_clear;
*push = *examplec::IntVector_push;
*pop = *examplec::IntVector_pop;
*get = *examplec::IntVector_get;
*set = *examplec::IntVector_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        examplec::delete_IntVector($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : example::DoubleVector ##############

package example::DoubleVector;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( example );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = examplec::new_DoubleVector(@_);
    bless $self, $pkg if defined($self);
}

*size = *examplec::DoubleVector_size;
*empty = *examplec::DoubleVector_empty;
*clear = *examplec::DoubleVector_clear;
*push = *examplec::DoubleVector_push;
*pop = *examplec::DoubleVector_pop;
*get = *examplec::DoubleVector_get;
*set = *examplec::DoubleVector_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        examplec::delete_DoubleVector($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package example;

*PGSd_EOS_AM = *examplec::PGSd_EOS_AM;
*PGSd_EOS_PM = *examplec::PGSd_EOS_PM;
1;
