package BSplines::Commons;

use strict;
use warnings;

use List::Util qw(sum);
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(build_knots findout_interval build_bases);


sub build_knots {
	my ($s,$m,$p) = @_;
	my $n = $#$s+1;
	my $u = [(0) x ($m+$p+1)];
	($$u[$_-$m],$$u[$_]) = ($$s[0],$$s[-1]) foreach ($m..$m+$p);
	$$u[$_] += (sum @{$s}[$_-$p-1..$_+$n-$m])/($n-$m+$p+2) foreach ($p+1..$m-1);
	return $u;
}

sub findout_interval {
	my ($u,$p,$up) = @_;
	my $i = $p;
	if ($up>=$$u[-1]) {
		$i = @$u-$p-2;
	} elsif ($up>$$u[0] and $up<$$u[-1]) {
		my $l = @$u;
		while ($l>1) {
			my $h = int($l/2);
			$i += $h if ($up>=$$u[$i+$h]);
			$l -= $h;
		}
	}
	return $i;
}

sub build_bases {
	my ($u,$p,$up) = @_;
	my $k = findout_interval $u,$p,$up;
	my $b = [(0) x ($p+1)];
	$$b[0] = 1.0;
	for (my $j=1; $j <= $p; $j++) {
		my $i = $k;
		$$b[$j] = $$b[$j-1]*($up-$$u[$i])/($$u[$i+$j]-$$u[$i]);
		for (my $l=$j-1; $l >= 1; $l--) {
			$i = $k+$l-$j;
			$$b[$l] = $$b[$l-1]*($up-$$u[$i])/($$u[$i+$j]-$$u[$i])+
          $$b[$l]*($$u[$i+$j+1]-$up)/($$u[$i+$j+1]-$$u[$i+1]);
		}
		$i = $k-$j;
		$$b[0] *= ($$u[$i+$j+1]-$up)/($$u[$i+$j+1]-$$u[$i+1]);
	}
	$k -= $p;
	return ($k,$b);
}

