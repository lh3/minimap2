#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>33088);
getopts('n:', \%opts);

my $pseudo = .5;
my $tot = $pseudo;
my $err = $pseudo;
while (<>) {
	chomp;
	if (/^Q\t(\d+)\t(\d+)\t(\d+)/) {
		$tot += $2;
		$err += $3;
		print join("\t", $1, $err/$tot, $tot / $opts{n}), "\n";
	}
}
