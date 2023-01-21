#!/usr/bin/env perl -w
# Simplest possible revcomper    
use strict;

my %substitute = (A => "T", T => "A", C => "G", G => "C", a => "T", t => "A", c => "G", g => "C");

my $seq = shift;
my $rev = "";
my $mend = "";
while ($seq ne "") {
    $mend = chop $seq;
    $rev .= $substitute{"$mend"} || $mend;
}
print "$rev";
