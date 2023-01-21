#!/usr/bin/env perl -w
# Wrangledups - 
# Usage -- perl wrangleDups.pl -r crclass -d 
#Created Tuesday, 29 November 2022.
#    Carrie Ganote 
    
use strict;
use Getopt::Long;
use Data::Dumper qw(Dumper);

my %substitute = (A => "T", T => "A", C => "G", G => "C", a => "T", t => "A", c => "G", g => "C");

my $output = "test.out";
my $crc = "test.crclass";
my $dup = "testduplicates";
# Threshold means a minimum spacer length
my $threshold = 16;

GetOptions("crc|r=s" => \$crc,
	   "dup|d=s" => \$dup,
	   "threshold|t=i" => \$threshold,
           "output|o=s" => \$output) or die ("No dice with options: $!");

sub revcomp{
    my $seq = shift;
    my $rev = "";
    my $mend = "";
    while ($seq ne "") {
        $mend = chop $seq;
	$rev .= $substitute{"$mend"} || $mend;
    }
    return "$rev";
}

my %drclass;
my @dups;
my $countDown = 30;

open (my $crcfh, $crc) or die "Error, cannot open file $crc";
my $header = <$crcfh>;
while (<$crcfh>) {
    chomp;
    my @hit = split(/,/);
    # An example line: ATCTTGATCACATCCTATCGAAATTTCC,VI-B,0.0271991677582264,Clostridium_botulinum_A634,NZ_CP013844.1,33103..33654,I-B,ATTTAAAAACATCATATGCTATATTCA,10
   # print "$hit[0] $hit[2] $hit[-1] $hit[1]==$hit[6]\n";
    if ($hit[2] > 0.75 || $hit[-1] <= 4 && $hit[1] eq $hit[6]){
#	print "\t match! $hit[0] $hit[2] $hit[-1]\n";
	$drclass{$hit[0]} = sprintf( "$hit[1] %.3f",  $hit[2]) ;
	$drclass{revcomp($hit[0])} = sprintf ("$hit[1] %.3f",  $hit[2]);
    }
}
#print Dumper(\%drclass);
close $crcfh;

open (my $dupfh, $dup) or die "Error, cannot open file $dup";
my $lastheader = "";
my $lastrepeat = "";
my $lastsample = "";
## These are hashes of arrays
my %spacers_in_samples;
my %spacers_in_repeats;
my %samples_in_repeats;

while (<$dupfh>) {
    chomp;
    # An example: >A7_02_G29DR1 \n CGGTTCATCCCCGCGTAGGCGGGGAACA	\n >A7_02_G29SP6_Cov_1
    if (/^>/){
	$lastheader = $_;
	($lastsample) = />(.._\d+)/;
    }
    else {
	if ($lastheader =~ m/_Cov_/){ # it is a spacer
	    if ($drclass{$lastrepeat} && length($_) > $threshold){
		push @{$spacers_in_samples{$_}}, $lastsample;
		push @{$spacers_in_samples{revcomp($_)}}, $lastsample;
		
		push @{$spacers_in_repeats{$_}}, $lastrepeat;
		push @{$spacers_in_repeats{revcomp($_)}}, $lastrepeat;
		#push @{$spacers_in_repeats{$_}}, revcomp($lastrepeat);
		#push @{$spacers_in_repeats{revcomp($_)}}, revcomp($lastrepeat);
	    }
	}
	elsif ($lastheader =~ m/DR/){ # it is a repeat
	    $lastrepeat = $_;
	    if ($drclass{$lastrepeat}){
		push @{$samples_in_repeats{$lastrepeat}}, $lastsample;
		push @{$samples_in_repeats{revcomp($lastrepeat)}}, $lastsample;
		push @{$samples_in_repeats{$lastrepeat}}, revcomp($lastsample);
		push @{$samples_in_repeats{revcomp($lastrepeat)}}, revcomp($lastsample);
	    }
	}
	else {
	    die("What is going on with this line? $lastheader $_\n");
	}
    }
}
#print "samples in repeats: \n";
#print Dumper(\%samples_in_repeats);
#print "spacers in samples: \n";
#print Dumper(\%spacers_in_samples);
#print "spacers in repeats: \n";
#print Dumper(\%spacers_in_repeats);
close $dupfh;
open (my $out, '>', $output) or print "Could not open file $output: $!";

my $id = 0;
lewp: foreach my $spacer (keys %spacers_in_samples){
    my $header = ">ID${id}_";
    $id++;
    $header .= join("__",  @{$spacers_in_samples{$spacer}});
    $header .= " with repeats ";    
    foreach my $repeat (@{$spacers_in_repeats{$spacer}}){
	$header .= " $repeat $drclass{$repeat} "
    }
    $header .= "\n";
    print $out "$header$spacer\n";
    my $count = scalar @{$spacers_in_samples{$spacer}};
    if ($count > 2){
	print "$header: $count samples had $spacer\n";
	print Dumper(\@{$spacers_in_samples{$spacer}});
	print "these had repeats:\n";
	foreach my $repeat (@{$spacers_in_repeats{$spacer}}){
	    print "repeat is $repeat $drclass{$repeat}\n";
	}
    }
}
close $out;

#print Dumper %unique;
