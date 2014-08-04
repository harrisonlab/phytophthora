#!/usr/bin/perl -w
use strict;
use CWD;

my $line;

for (@_) {
	my @nameparts = split ('/', $_);
	$genome = $nameparts[-3];
	$isolate = $nameparts[-2];
	open (INFILE, "$_") or die;
	while ($line = <INFILE>) {
		$line =~ s/>/>$genome\|$isolate\|/g;
		$line =~ s/\s//g;
		print $line;
	}
}

exit
		