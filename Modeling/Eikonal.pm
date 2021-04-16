package Modeling::Eikonal;

use strict;
use warnings;

use BSplines::Surface;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw();


sub new {
	my ($class,$model) = @_;
	return bless $model,$class;
}

sub initialize {
	my ($self,$source) = @_;
	$$self{time} = [map { [(999) x $$self{shape}[1]] } 1..$$self{shape}[0]];
	$$self{tag} = [map { [(-1) x $$self{shape}[1]] } 1..$$self{shape}[0]];
	my $ix = int (($$source[0]-$$self{bounds}[0])/$$self{increments}[0]);
	my $iz = int (($$source[1]-$$self{bounds}[2])/$$self{increments}[1]);
	if ($$source[0]==$$self{bounds}[0]+$ix*$$self{increments}[0] or 
			$$source[1]==$$self{bounds}[2]+$iz*$$self{increments}[1]) {
		$$self{time}[$ix][$iz] = 0.0;
		$$self{tag}[$ix][$iz] = 0;
		$$self{front} = [[$ix,$iz]];
	} else {
		my $model = {
				shape => [int($$self{increments}[0]/10.0)+1,int($$self{increments}[1]/10.0)+1],
				bounds => [$$source[0]-$$self{increments}[0],$$source[0]+$$self{increments}[0],
						$$source[1]-$$self{increments}[1],$$source[1]+$$self{increments}[1]],
			};
		$$model{increments} = [$$self{increments}[0]/($$model{shape}[0]-1),
				$$self{increments}[1]/($$model{shape}[1]-1)];
		my $x = [map { $$self{bounds}[0]+$$self{increments}[0]*$_ } 0..$$self{shape}[0]-1];
		my $z = [map { $$self{bounds}[2]+$$self{increments}[1]*$_ } 0..$$self{shape}[1]-1];
		my $v = new BSplines::Surface($x,$z,$$self{speed},$$self{shape}[0],$$self{shape}[1],3,3);
		for (my $ixx=0; $ixx < $$model{shape}[0]; $ixx++) {
			$$model{speed}[$ixx][$_] = $v->evaluate($$model{bounds}[0]+$$model{increments}[0]*$ixx,
			$$model{bounds}[2]+$$model{increments}[1]*$_) foreach (0..$$model{shape}[1]-1);
		}
		$$model{time} = [map { [(999) x $$model{shape}[1]] } 1..$$model{shape}[0]];
		$$model{tag} = [map { [(-1) x $$model{shape}[1]] } 1..$$model{shape}[0]];
		my $ixx = int (($$source[0]-$$model{bounds}[0])/$$model{increments}[0]);
		my $izz = int (($$source[1]-$$model{bounds}[2])/$$model{increments}[1]);
		$$model{time}[$ixx][$izz] = 0.0;
		$$model{tag}[$ixx][$izz] = 0;
		$$model{front} = [[$ixx,$izz]];
		$model = bless $model,__PACKAGE__;
		$model->move_node($model->remove_node()) while (@{$$model{front}}>0);
		my $xx = [map { $$model{bounds}[0]+$$model{increments}[0]*$_ } 0..$$model{shape}[0]-1];
		my $zz = [map { $$model{bounds}[2]+$$model{increments}[1]*$_ } 0..$$model{shape}[1]-1];
		my $t = new BSplines::Surface($xx,$zz,$$model{time},$$model{shape}[0],$$model{shape}[1],3,3);
		$$self{time}[$ix][$iz] = $t->evaluate($$x[$ix],$$z[$iz]);
		$$self{tag}[$ix][$iz] = 0;
		$$self{time}[$ix+1][$iz] = $t->evaluate($$x[$ix+1],$$z[$iz]);
		$$self{tag}[$ix+1][$iz] = 0;
		$$self{time}[$ix][$iz+1] = $t->evaluate($$x[$ix],$$z[$iz+1]);
		$$self{tag}[$ix][$iz+1] = 0;
		$$self{time}[$ix+1][$iz+1] = $t->evaluate($$x[$ix+1],$$z[$iz+1]);
		$$self{tag}[$ix+1][$iz+1] = 0;
		$$self{front} = [[$ix,$iz],[$ix+1,$iz],[$ix,$iz+1],[$ix+1,$iz+1]];
	}
}

sub remove_node {
	my ($self) = @_;
	$$self{front} = [reverse sort { $$self{time}[$$a[0]][$$a[1]]<=>
			$$self{time}[$$b[0]][$$b[1]] } @{$$self{front}}];
	my $node = pop @{$$self{front}};
	$$self{tag}[$$node[0]][$$node[1]] = 1;
	return $node;
}

sub insert_node {
	my ($self,$node) = @_;
	$$self{tag}[$$node[0]][$$node[1]] = 0;
	push @{$$self{front}},$node;
}

sub has_node {
	my ($self,$node) = @_;
	return ($$node[0]>=0 and $$node[0]<$$self{shape}[0]
			and $$node[1]>=0 and $$node[1]<$$self{shape}[1]);
}

sub get_neighbors {
	my ($self,$node) = @_;
	my @neighbors = (
			[$$node[0]-1,$$node[1]],
			[$$node[0]+1,$$node[1]],
			[$$node[0],$$node[1]-1],
			[$$node[0],$$node[1]+1]
		);
	return grep { $self->has_node($_) } @neighbors;
}

sub compute_time {
	my ($self,$node) = @_;
	my @neighbors = grep { $$self{tag}[$$_[0]][$$_[1]]==1
			} $self->get_neighbors($node);
	my ($a,$b,$c) = (0,0,-$$self{speed}[$$node[0]][$$node[1]]**(-2));
	my $kappa = [(0) x 3];
	foreach my $node_a (@neighbors) {
		my ($h,$node_b);
		if ($$node[0] == $$node_a[0]) {
			$h = $$self{increments}[1];
			$node_b = [$$node[0],2*$$node_a[1]-$$node[1]];
		} elsif ($$node[1] == $$node_a[1]) {
			$h = $$self{increments}[0];
			$node_b = [2*$$node_a[0]-$$node[0],$$node[1]];
		}
		if ($self->has_node($node_b) and
				$$self{tag}[$$node_b[0]][$$node_b[1]]==1) {
			$$kappa[0] = 2.250*$h**(-2);
			$$kappa[1] = 0.666667*$$kappa[0]*($$self{time}[$$node_b[0]][$$node_b[1]]-
					4.0*$$self{time}[$$node_a[0]][$$node_a[1]]);
			$$kappa[2] = 0.111111*($h*$$kappa[1])**2;
		} else {
			$$kappa[0] = $h**(-2);
			$$kappa[1] = -2.0*$$kappa[0]*$$self{time}[$$node_a[0]][$$node_a[1]];
			$$kappa[2] = 0.25*$$kappa[1]**2/$$kappa[0];
		}
		$a += $$kappa[0];
		$b += $$kappa[1];
		$c += $$kappa[2];
	}
	my $delta = abs($b**2-4*$a*$c);
	$$self{time}[$$node[0]][$$node[1]] = 0.5*(-$b+sqrt($delta))/$a;
}

sub move_node {
	my ($self,$node) = @_;
	my @neighbors = $self->get_neighbors($node);
	$self->insert_node($_) foreach (grep { $$self{tag}[$$_[0]][$$_[1]]==-1 } @neighbors);
	$self->compute_time($_) foreach (grep { $$self{tag}[$$_[0]][$$_[1]]==0 } @neighbors);
}

sub solve {
	my ($self,$source) = @_;
	$self->initialize($source);
	$self->move_node($self->remove_node()) while (@{$$self{front}}>0);
}

1;
