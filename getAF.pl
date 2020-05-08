#!/usr/bin/perl -w
use strict;

while(<>){
  next if (/^#/);
  my @f = split(/\s+/);
  my $xtra = $f[9];
  my @g = split(/:/, $xtra);
  my $r = $g[2];
  print "$r\n" if ($r);
}
