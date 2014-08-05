#!/usr/bin/perl -w
use strict;
use Cwd;

my $usage = "append_rxlr.pl <infile_sp_rxlr.fa> > appended_rxlr_sp.fa\n";
my $line;

for (@ARGV) {
	my @nameparts = split ('/', $_);
	my $genome = $nameparts[-3];
	my $isolate = $nameparts[-2];
	my $file = $_;
	open (INFILE, "$file") or die $usage;
	while ($line = <INFILE>) {
		$line =~ s/>/>$genome\|$isolate\|/g;
		$line =~ s/\s//g;
		print "$line\n";
	}
	close (INFILE);
}

exit
		