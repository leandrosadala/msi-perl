=pod

=head1

=head1 NAME

MSI-INPUT-CREATION

=head1 DESCRIPTION

Creates json input file for msi application

=head1 SYNOPSIS

msi-input-creation -m MODEL -b BOUNDS -i INCREMENTS -o OUTPUT 

                   [-l LEVEL] [-j JOBS]

=head1 REQUIREMENTS

=over 4

=item B<-m,--model>

json velocity model file

=item B<-b,--bounds>

bounds of the image

=item B<-i,--increments>

increments of the image

=item B<-o,--output>

json created file

=back

=head1 OPTIONS

=over 4

=item B<-l,--level>

registration level of the records

=item B<-j,--jobs>

no of concurrent processes

=back

=head1 EXAMPLES

=over 4

=item registration level set to the lower depth model bound

msi-input-creation -m model.json -b 0,10000,0,2000 -i 10,10

=item registration level set to a given depth

msi-input-creation -m model.json -b 0,10000,0,2000 -i 10,10 -l 10

=item registration level set to a given depth with 4 jobs

msi-input-creation -m model.json -b 0,10000,0,2000 -i 10,10 -l 10 -j 4

=back

=head1 CREDITS

Leandro Sadala,L<<leandrosadala@petrobras.com.br>>

=cut

use strict;
use warnings;

use JSON;
use File::Basename;
use Modeling::FirstArrival;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use MCE::Map;
use Data::Dump qw(pp);

sub construct_curves {
	my ($index,$image,$model,$level) = @_;
	
    my $ix = int($index/$$image{shape}[1]);
    my $iz = $index % $$image{shape}[1];
    
    my $source = [$$image{bounds}[0]+$$image{increments}[0]*$ix,
            $$image{bounds}[2]+$$image{increments}[1]*$iz];
    
    my $es = new Modeling::Eikonal($model);
    $es->solve($source);
    
    $$model{t} = new BSplines::Surface($$model{x},$$model{z},$$model{time},$$model{shape}[0],$$model{shape}[1],3,3);
    my $fstarr = new Modeling::FirstArrival($$model{v},$$model{t});
    
    my $receivers = [map { [$_,$level] } @{$$model{x}}];
    
    my @curves = $fstarr->build($source,$receivers);
    
    print ">> first arrival of source ($$source[0] m,$$source[1] m) built <<\n";

	return \@curves;
	
}

sub main {
  
  my $options = {
      'help' => undef,
      'model' => undef,
      'bounds' => undef,
      'output' => undef,
      'increments' => undef,
      'level' => undef,
      'jobs' => 1,
    };
  GetOptions($options,'help|h!','model|m=s','bounds|b=s',
      'output|o=s','increments|i=s','level|l=f','jobs|j=i');
  pod2usage(-exitval => 0, -verbose => 2) if $$options{help};
  pod2usage(-exitval => 1, -verbose => 1) unless ($$options{model} 
      and $$options{bounds} and $$options{output} or $$options{increments});
  
  open IN,"<$$options{model}";
  read IN,my $expr,-s IN;
  close IN;
  
  my $model = JSON->new->decode($expr);
  
  $$options{level} = $$model{bounds}[2] unless $$options{level};
  
  open IN,'<',join('/',dirname($$options{model}),$$model{filepath});
  foreach (0..$$model{shape}[0]-1) {
    read IN,$expr,4*$$model{shape}[1];
    $$model{speed}[$_] = [unpack "f$$model{shape}[1]",$expr];
  }
  close IN;
  
  $$model{x} = [map { $$model{increments}[0]*$_+$$model{bounds}[0] } 0..$$model{shape}[0]-1];
  $$model{z} = [map { $$model{increments}[1]*$_+$$model{bounds}[2] } 0..$$model{shape}[1]-1];
  $$model{v} = new BSplines::Surface($$model{x},$$model{z},$$model{speed},$$model{shape}[0],$$model{shape}[1],3,3);
  
  my $image = {
      bounds => [map { 1.0*$_ } split(/,/,$$options{bounds})],
      increments => [map { 1.0*$_ } split(/,/,$$options{increments})],
      filepath => basename($$options{output}),
    };
  
  $$image{'first arrivals'}{positions} = $$model{x};
  
  $$image{shape} = [int(($$image{bounds}[1]-$$image{bounds}[0])/$$image{increments}[0])+1,
      int(($$image{bounds}[3]-$$image{bounds}[2])/$$image{increments}[1])+1];
  
  MCE::Map::init { max_workers => $$options{jobs}, chunk_size => 1 };
  
  my @curves = mce_map { construct_curves $_,$image,$model,
        $$options{level} } 0..$$image{shape}[0]*$$image{shape}[1]-1;
  
  MCE::Map::finish;
  
  foreach (0..$#curves) {
  	$$image{'first arrivals'}{traveltimes}[$_] = $curves[$_][0];
  	$$image{'first arrivals'}{amplitudes}[$_] = $curves[$_][1];
  }

  open OUT,">$$options{output}";
  print OUT JSON->new->pretty(1)->encode($image);
  close OUT;

}

main;

