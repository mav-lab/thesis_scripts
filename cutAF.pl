#!/usr/bin/perl -w
use strict;
my $cutoff = shift || 0.3;

while(<>){
  if (/^#/){
    if (/^#[^#]/){
      print "##cutAF $cutoff\n";
    }
    print;
    next;
  }
  my @f = split(/\s+/);
  my $xtra = $f[9];
  my @g = split(/:/, $xtra);
  my $r = $g[2];
  if ($r < $cutoff){
    print;
  }
}
