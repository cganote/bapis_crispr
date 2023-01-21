#!/usr/bin/env perl -w
# scrape repeats from raw reads, given matching reads as input 
# Usage -- perl scrapeBombellaRepeats.pl -r "ATCG" -i test.txt
#Created Tuesday, 29 November 2022.
#    Carrie Ganote 
    
use strict;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use String::Approx 'adist';

my $input="STDIN";
my $output="STDOUT";

my %substitute = (A => "T", T => "A", C => "G", G => "C", a => "T", t => "A", c => "G", g => "C");

my $mismatches = 0;
my $minimum = 2;
my $repeat = "ATCG";

GetOptions("input|i=s" => \$input,
	   "repeat|r=s" => \$repeat,
	   "mismatches|n=i" => \$mismatches,
	   "minimum|m=s" => \$minimum,
           "output|o=s" => \$output) or die ("No dice with options: $!");

my $infh = *STDIN;
if ($input ne "STDIN"){
    open ($infh, $input) or die "Error, cannot open file $input";
}
my @mids;
my @lefts;
my @rights;
while (<$infh>) {
    chomp;
    my @chunks = split(/$repeat/, $_, -1);
    push (@rights, pop (@chunks));
    push (@lefts, shift(@chunks));
    push (@mids, @chunks);
}
close $infh;

# get size of mids. Remember - the minimum mid is a good starting point. There will be a breakpoint
# where this size + the size of the repeat indicates a repeat error and at least two mids
my $midsize = 100;
my $replen = length($repeat);
foreach my $mid (@mids){
    my $len = length($mid);
    $midsize = $len if $len < $midsize;
}
sub chewBackward{
    my @splitRA = @_;
    foreach my $go (@splitRA){
	## If there is no repeat in line, the left (backward) array will be undef. Check for $go first
	if ($go && length($go) > $minimum){
	    my $pretty = "";
	    while (length($go) > $midsize){
		$pretty = substr($go,-$midsize) . " " . $pretty;
		$go = substr($go,0,-$midsize);
		#print "go is now $go\n";
		if (length($go) > $replen){
		    $pretty = substr($go,-$replen) . " " . $pretty;
		    $go = substr($go,0,-$replen);
		    #print "   go is now $go\n";
		}
	    }
	    print "$go $pretty\n";
	}
    }
}
sub chewForward{
    my @splitRA = @_;
    foreach my $go (@splitRA){
	if (length($go) > $minimum){
	    my $pretty = "";
	    while (length($go) > $midsize){
		$pretty .= substr($go,0,$midsize) . " "; 
		$go = substr($go,$midsize);	    
		#print "go is now $go\n";
		if (length($go) > $replen){
		    $pretty .= substr($go,0,$replen) . " ";
		    $go = substr($go,$replen);		
		    #print "    go is now $go\n";
		}
	    }
	    print "$pretty$go\n";
	}
    }
}
print "Rights:\n";
#print Dumper(\@rights);
chewForward(@rights);
print "Mids:\n";
#print Dumper(\@mids);
chewForward(@mids);
print "Lefts:\n";
#print Dumper(\@lefts);
chewBackward(@lefts);

## Todo: Dedup, find perfect substrings and collapse, 
#    extend (seed by minimum on each side) to create consensus.

## This is experimental. I'm gonna ignore it for now for time reasons,
#   and just use split and chew and clean up by hand
sub seek{
    my ($string, $repeat, $size) = @_;    
    my $i = 0;
    my $j = 0;
    my $maxscore = 0;
    my $misses = $mismatches;
    # Start at first character of our repeat. Life is easy if it matches
    while ($j < length($repeat)){
	my $complete = "";
	my $score = 0;
	my $startingI = $i;
	# I feel like we could start up to $mismatches away from where we expect the next rep to be
	while ($i < ($repeat + $size)){
	    if (substr($string, $i, 1) eq substr($repeat, $j, 1)){
		$score++;
	    }
	    else{
		$mismatches++;
	    }
	    $i++;
	}
	$j++;
    }
}

