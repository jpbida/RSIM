#!/usr/bin/perl
use Math::VectorReal;
use Math::Trig;
use constant PI    => 4 * atan2(1, 1);
@p1=(0,0,0);  
@p2=(1,0,0); 
@p3=(1,1,0); 
@p4=(0,1,0); 

my $hash=();
$ang=0;
while($ang < 2*PI){
$hash->{'x'}=$p4[0];
$hash->{'y'}=$p4[1];
$hash->{'z'}=$p4[2];
$hash=translate($hash,1,0,0);
$hash=rotateY($hash,$ang);
$hash=translate($hash,-1,0,0);
my $p4i=();
$p4i[0] =$hash->{'x'};
$p4i[1] =$hash->{'y'};
$p4i[2] =$hash->{'z'};
print "$p4i[0],$p4i[1],$p4i[2]\n";
$ang=$ang+PI/2;
tor(\@p1,\@p2,\@p3,\@p4i);

}

sub tor{
#### Calculate Torsion Angles ####
use Math::Trig;
use Math::VectorReal;
my	$p1=shift(@_);
my	$p2=shift(@_);
my	$p3=shift(@_);
my	$p4=shift(@_);
my @v1=@$p1;
my  @v2=@$p2;
my  @v3=@$p3;
my  @v4=@$p4;

#### Create vectors ####
my $a=vector(($v2[0]-$v1[0]),($v2[1]-$v1[1]),($v2[2]-$v1[2]));
my $b=vector(($v3[0]-$v2[0]),($v3[1]-$v2[1]),($v3[3]-$v2[2]));
my $c=vector(($v3[0]-$v4[0]),($v3[1]-$v4[1]),($v3[2]-$v4[2]));

### Normal vectors to the two planes ###
my $n1=$a x $b;
my $n2=$b x $c;
my $tor = acos(($n1.$n2)/($n1->length * $n2->length));
#### Determining the sign of the angle ####
my $sv = acos($n1.$c/($n1->length*$c->length));
$m=PI;
$m=$m/2;
my $sign=0;
if($sv <= $m){
$sign = 1;}else{
$sign = -1;
}
$tor = $tor*$sign;
print "Tor:$tor\n";


}

sub translate{
my $hash=shift(@_);
my $tx=shift(@_);
my $ty=shift(@_);
my $tz=shift(@_);
$hash->{"x"}=$hash->{"x"}-$tx;
$hash->{"y"}=$hash->{"y"}-$ty;
$hash->{"z"}=$hash->{"z"}-$tz;
return $hash;
}

sub rotateY{
my $hash=shift(@_);
my $ang=shift(@_);
$z1=$hash->{"z"};
$x1=$hash->{"x"};
#	about Y
$hash->{"x"}=$x1*cos($ang)+$z1*sin($ang);
$hash->{"z"}=$z1*cos($ang)-$x1*sin($ang);
return $hash;
}
sub rotateZ{
my $hash=shift(@_);
my $ang=shift(@_);
$y1=$hash->{"y"};
$x1=$hash->{"x"};
#	about Y
$hash->{"x"}=$x1*cos($ang)+$y1*sin($ang);
$hash->{"y"}= -1*$x1*sin($ang)+$y1*cos($ang);
return $hash;
}

sub rotateX{
	use Math::Trig;
my $hash=shift(@_);
my $ang=shift(@_);
$z1=$hash->{"z"};
$y1=$hash->{"y"};
#	about 

$hash->{"y"}=$y1*cos($ang)-$z1*sin($ang);
$hash->{"z"}=$y1*sin($ang)+$z1*cos($ang);
return $hash;
}

