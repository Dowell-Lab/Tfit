#!/usr/bin/perl

my $cnt = 0;
while(<>) {
  if (! /^#/) {$cnt++; }
}
print $cnt ."\n";
