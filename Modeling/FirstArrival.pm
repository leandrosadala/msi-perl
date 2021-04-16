package Modeling::FirstArrival;

use strict;
use warnings;

use BSplines::Surface;
use Modeling::Eikonal;
use MCE::Map;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();


sub new {
	my ($class,@args) = @_;
	return bless {v=>$args[0],t=>$args[1]},$class;
}

sub step_backward {

	my ($self,$r,$h) = @_;

	my ($k1r,$k2r,$k3r,$k4r);
	my ($vr,$gr,$pr);

	my $rp = [@$r];

	$vr = $$self{v}->evaluate($$r[0],$$r[1]);
	$gr = $$self{t}->gradient($$r[0],$$r[1]);
	$pr = [-$$gr[0],-$$gr[1]];
	$k1r = [$$pr[0]*$vr**2,$$pr[1]*$vr**2];

	$r = [$$rp[0]+0.5*$h*$$k1r[0],$$rp[1]+0.5*$h*$$k1r[1]];
	$vr = $$self{v}->evaluate($$r[0],$$r[1]);
	$gr = $$self{t}->gradient($$r[0],$$r[1]);
	$pr = [-$$gr[0],-$$gr[1]];
	$k2r = [$$pr[0]*$vr**2,$$pr[1]*$vr**2];
	
	$r = [$$rp[0]+0.5*$h*$$k2r[0],$$rp[1]+0.5*$h*$$k2r[1]];
	$vr = $$self{v}->evaluate($$r[0],$$r[1]);
	$gr = $$self{t}->gradient($$r[0],$$r[1]);
	$pr = [-$$gr[0],-$$gr[1]];
	$k3r = [$$pr[0]*$vr**2,$$pr[1]*$vr**2];
	
	$r = [$$rp[0]+$h*$$k3r[0],$$rp[1]+$h*$$k3r[1]];
	$vr = $$self{v}->evaluate($$r[0],$$r[1]);
	$gr = $$self{t}->gradient($$r[0],$$r[1]);
	$pr = [-$$gr[0],-$$gr[1]];
	$k4r = [$$pr[0]*$vr**2,$$pr[1]*$vr**2];
	
	$r = [
	    $$rp[0]+$h*($$k1r[0]+2.0*($$k2r[0]+$$k3r[0])+$$k4r[0])/6.0,
			$$rp[1]+$h*($$k1r[1]+2.0*($$k2r[1]+$$k3r[1])+$$k4r[1])/6.0
		];
	
	return $r;
	
}

sub compute_amplitude {
	my ($self,$r,$tau) = @_;
	my ($h,$l,$dl,$rp) = (5e-3,0.0);
	while ($tau>0.0) {
		$rp = [@$r];
		$r = $self->step_backward($r,$h);
		$dl = sqrt(($$r[0]-$$rp[0])**2+($$r[1]-$$rp[1])**2);
		$l += $dl;
		$tau -= $h;
	}
	if ($tau<0.0) {
		$l -= $dl;
		my $vr = $$self{v}->evaluate($$rp[0],$$rp[1]);
		$$r[0] = $$rp[0]+($$r[0]-$$rp[0])*($tau+$h)/$h;
		$$r[1] = $$rp[1]+($$r[1]-$$rp[1])*($tau+$h)/$h;
		$dl += sqrt(($$r[0]-$$rp[0])**2+($$r[1]-$$rp[1])**2);
		$l += $dl;
		$tau += $dl/$vr-$h;
	}
	return 1e3/($l*sqrt($l));
}

sub build {
	my ($self,$s,$r) = @_;
	my $time = [map { $$self{t}->evaluate($$_[0],$$_[1]) } @$r];
	my $ampl = [mce_map { $self->compute_amplitude($$r[$_],$$time[$_]) } 0..$#$r];
	return ($time,$ampl);
}

1;
