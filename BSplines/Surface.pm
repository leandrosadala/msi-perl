package BSplines::Surface;

use strict;
use warnings;

use Math::Matrix;
use BSplines::Commons;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();


sub build_control_points {
	my ($s,$t,$D,$u,$v,$m,$n,$p,$q) = @_;
	my $Q = new Math::Matrix($D)->swaprc;
	my $N = Math::Matrix->zeros($#$t+1,$n);
	for (my $i=0; $i < @$t; $i++) {
		my ($k,$b) = build_bases $v,$q,$$t[$i];
		$$N[$i][$k+$_] = $$b[$_] foreach (0..$#$b);
	}
	my $Nt = $N->swaprc;
	$Q = $Nt->mmul($Q);
	my $A = $Nt->mmul($N);
	$Q = $A->concat($Q)->solve;
	$Q = $Q->swaprc;
	$N = Math::Matrix->zeros($#$s+1,$m);
	for (my $i=0; $i < @$s; $i++) {
		my ($k,$b) = build_bases $u,$p,$$s[$i];
		$$N[$i][$k+$_] = $$b[$_] foreach (0..$#$b);
	}
	$Nt = $N->swaprc;
	my $P = $Nt->mmul($Q);
	$A = $Nt->mmul($N);
	$P = $A->concat($P)->solve;
	return $P->as_array;
}

sub evaluate {
	my ($self,$up,$vp) = @_;
	my ($k,$a) = build_bases $$self{u},$$self{p},$up;
	my ($l,$b) = build_bases $$self{v},$$self{q},$vp;
	my $value = 0.0;
  for (my $i=0; $i<=$#$a; $i++) {
    $value += $$a[$i]*$$b[$_]*$$self{P}[$k+$i][$l+$_] foreach (0..$#$b);
  }
	return $value;
}

sub gradient {
	my ($self,$up,$vp) = @_;
	my $gradient;
	my ($k,$a) = build_bases [@{$$self{u}}[1..$#{$$self{u}}-1]],$$self{p}-1,$up;
	my ($l,$b) = build_bases $$self{v},$$self{q},$vp;
	$$gradient[0] = 0.0;
	for (my $i=0; $i <=$#$a ; $i++) {
	  $$gradient[0] += $$a[$i]*$$b[$_]*($$self{P}[$k+$i+1][$l+$_]-
        $$self{P}[$k+$i][$l+$_])/($$self{u}[$k+$i+$$self{p}+1]-
          $$self{u}[$k+$i+1]) foreach (0..$#$b);
	}
	$$gradient[0] *= $$self{p};
	($k,$a) = build_bases $$self{u},$$self{p},$up;
	($l,$b) = build_bases [@{$$self{v}}[1..$#{$$self{v}}-1]],$$self{q}-1,$vp;
	$$gradient[1] = 0.0;
	for (my $i=0; $i <=$#$a ; $i++) {
	  $$gradient[1] += $$a[$i]*$$b[$_]*($$self{P}[$k+$i][$l+$_+1]-
          $$self{P}[$k+$i][$l+$_])/($$self{v}[$l+$_+$$self{q}+1]-
            $$self{v}[$l+$_+1]) foreach (0..$#$b);
	}
	$$gradient[1] *= $$self{q};
	return $gradient;
}

sub new {
	my ($class,@args) = @_;
	my $self = { 'm' => $args[3], 'n' => $args[4], 
	    'p' => $args[5], 'q' => $args[6] };
	$$self{u} = build_knots $args[0],$args[3],$args[5];
	$$self{v} = build_knots $args[1],$args[4],$args[6];
	$$self{P} = build_control_points $args[0],$args[1],$args[2],
			$$self{u},$$self{v},$args[3],$args[4],$args[5],$args[6];
	return bless $self,$class;
}

1;
