#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>33088, s=>100);
getopts('n:', \%opts);

my $pseudo = .5;
my $tot = $pseudo;
my $err = $pseudo;
my $tot_last_out = -$opts{s};
my $state = 0;
my $mapq = 0;
while (<>) {
	chomp;
	if (/^Q\t(\d+)\t(\d+)\t(\d+)/) {
		$tot += $2;
		$err += $3;
		if ($tot - $tot_last_out >= $opts{s}) {
			print join("\t", $1, $err/$tot, $tot / $opts{n}), "\n";
			$tot_last_out = $tot;
			$state = 0;
		} else {
			$state = 1;
			$mapq = $1;
		}
	}
}
if ($state) {
	print join("\t", $mapq, $err/$tot, $tot / $opts{n}), "\n";
}
