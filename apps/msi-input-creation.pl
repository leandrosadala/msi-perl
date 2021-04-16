=pod

=head1

=head1 NAME

MSI-INPUT-CREATION

=head1 DESCRIPTION

Creates json input file for msi application

=head1 SYNOPSIS

msi-input-creation -m=MODEL -b=BOUNDS -i=INCREMENTS -o=OUTPUT [-l=LEVEL]

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

=item B<-l,-level>

registration level of the records

=back

=head1 EXAMPLES

=over 4

=item registration level set to the lower depth model bound

msi-input-creation -m=model.json -b=0,10000,0,2000 -i=10,10

=item registration level set to a given depth

msi-input-creation -m=model.json -b=0,10000,0,2000 -i=10,10 -l=10

=back

=head1 CREDITS

Leandro Sadala,L<leandrosadala@petrobras.com.br>

=cut

use strict;
use warnings;

use JSON;
use File::Basename;
use Modeling::FirstArrival;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);

sub main {
  
  my $options = {
      'help' => undef,
      'model' => undef,
      'bounds' => undef,
      'output' => undef,
      'increments' => undef,
      'level' => undef,
    };
  GetOptions($options,'help|h!','model|m=s','bounds|b=s',
      'output|o=s','increments|i=s','level|l=f');
  pod2usage(-exitval => 0, -verbose => 2) if $$options{help};
  pod2usage(-exitval => 1, -verbose => 1) unless ($$options{model} 
      and $$options{bounds} and $$options{output} or $$options{increments});
  
  open IN,"<$$options{model}";
  read IN,my $expr,-s IN;
  close IN;
  
  my $model = JSON->new->decode($expr);
  
  $$options{level} = $$model{bounds}[2] unless $$options{level};
  
  my $v;
  open IN,'<',join('/',dirname($$options{model}),$$model{filepath});
  foreach (0..$$model{shape}[0]-1) {
    read IN,$expr,4*$$model{shape}[1];
    $$v[$_] = [unpack "f$$model{shape}[1]",$expr];
  }
  close IN;
  
  $$model{speed} = $v;
  
  my $x = [map { $$model{increments}[0]*$_+$$model{bounds}[0] } 0..$$model{shape}[0]-1];
  my $z = [map { $$model{increments}[1]*$_+$$model{bounds}[2] } 0..$$model{shape}[1]-1];
  
  $v = new BSplines::Surface($x,$z,$v,$$model{shape}[0],$$model{shape}[1],3,3);
  
  my $image = {
      bounds => [map { 1.0*$_ } split(/,/,$$options{bounds})],
      increments => [map { 1.0*$_ } split(/,/,$$options{increments})],
      filepath => basename($$options{output}),
    };
  
  $$image{'first arrivals'}{positions} = $x;
  
  $$image{shape} = [int(($$image{bounds}[1]-$$image{bounds}[0])/$$image{increments}[0])+1,
      int(($$image{bounds}[3]-$$image{bounds}[2])/$$image{increments}[1])+1];
  
  my $eikonal = new Modeling::Eikonal($model);

  foreach my $ix (0..$$image{shape}[0]-1) {
    my $source;
    $$source[0] = $$image{bounds}[0]+$$image{increments}[0]*$ix;
    foreach my $iz (0..$$image{shape}[1]-1) {
      my $index = $ix*$$image{shape}[1]+$iz;
      $$source[1] = $$image{bounds}[2]+$$image{increments}[1]*$iz;
      $eikonal->solve($source);
      my $t = $$model{time};
      $t = new BSplines::Surface($x,$z,$t,$$model{shape}[0],$$model{shape}[1],3,3);
      my $fa = new Modeling::FirstArrival($v,$t);
      my $receivers = [map { [$_,$$options{level}] } @$x];
      ($$image{'first arrivals'}{traveltimes}[$index],
          $$image{'first arrivals'}{amplitudes}[$index]) = 
              $fa->build($source,$receivers);
      print ">> first arrival of source ($$source[0] m,$$source[1] m) built <<\n";
    }
  }
  
  open OUT,">$$options{output}";
  print OUT JSON->new->pretty(1)->encode($image);
  close OUT;

}

main;

