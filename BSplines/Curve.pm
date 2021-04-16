package BSplines::Curve;

use strict;
use warnings;


use Math::Matrix;
use BSplines::Commons;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();


sub build_control_points {
	my ($s,$D,$u,$m,$p) = @_;
	my $P = new Math::Matrix($D)->swaprc;
	my $N = Math::Matrix->zeros($#$s+1,$m);
	for (my $i=0; $i < @$s; $i++) {
		my ($k,$b) = build_bases $u,$p,$$s[$i];
		$$N[$i][$k+$_] = $$b[$_] foreach (0..$#$b);
	}
	my $Nt = $N->swaprc;
	$P = $Nt->mmul($P);
	my $A = $Nt->mmul($N);
	$P = $A->concat($P)->solve->swaprc;
	return $$P[0];
}

sub evaluate {
	my ($self,$up) = @_;
	my ($k,$b) = build_bases $$self{u},$$self{p},$up;
	my $value = 0.0;
	$value += $$b[$_]*$$self{P}[$k+$_] foreach (0..$#$b);
	return $value;
}

sub derivate {
	my ($self,$up) = @_;
	my ($k,$b) = build_bases 
			[@{$$self{u}}[1..$#{$$self{u}}-1]],$$self{p}-1,$up;
	my $value = 0.0;
	$value += $$b[$_]*($$self{P}[$k+$_+1]-
			$$self{P}[$k+$_])/($$self{u}[$k+$_+$$self{p}+1]-
				$$self{u}[$k+$_+1]) foreach (0..$#$b);
	$value *= $$self{p};
	return $value;
}

sub new {
	my ($class,@args) = @_;
	my $self = { 'm' => $args[2], 'p' => $args[3] };
	$$self{u} = build_knots $args[0],$args[2],$args[3];
	$$self{P} = build_control_points $args[0],$args[1],
			$$self{u},$$self{m},$$self{p};
	return bless $self,$class;
}

1;
