#!/usr/bin/perl

#----------------------------
# Script to extract genes with hmmer hits
#----------------------------


use strict;
use warnings;
my $usage = "filter_gff_StartStop.pl <gene_models.gff>";
my $hmmStart = "------- ------ -----    ------- ------ -----   ---- --  --------";
my $printLine = 0;
my @hmmLines;
my $header;
my %fastaHash;

my $hmmFile = shift or die $usage;
my $fastaFile = shift or die $usage;


# open hmm file
open HMM, $hmmFile;

# collect hmm hits into array
# Only take the lines from the file that show the hit stats
while (my $line = <HMM>) {
	if ($line =~ m/^\s$/ and $printLine == 1) {last;}
	elsif ($line =~ m/------ inclusion threshold ------/) {next;}
	elsif ($printLine == 1) {push @hmmLines, $line }
	elsif ($line =~ m/$hmmStart/) {$printLine = 1;}
}

# collect fasta hits into dictionary
open FASTA, $fastaFile;
while (my $line = <FASTA>) {
	$line =~ s/^\s+|\s+$//g;
	if ($line =~ m/^>/) {$header = substr($line, 1)}
	else { $fastaHash{$header} .= $line; }
}

foreach (@hmmLines) {
	my $hmmHit = $_;
	$hmmHit =~ s/^\s+|\s+$//g;
	my @aoHit = split(/\s+/, $hmmHit);
	my $eval = shift @aoHit; 
	my $score = shift @aoHit; 
	my $bias = shift @aoHit;  
	my $best_eval = shift @aoHit; 
	my $best_score = shift @aoHit; 
	my $best_bias = shift @aoHit; 
	my $exp_domains = shift @aoHit; 
	my $num_domains = shift @aoHit; 
	my $seq_name = shift @aoHit; 
	my $description = shift @aoHit; 
# search for hmm hit names in dictionary
# if hit print hmm hit name, score, number of domains and e-number
# print hit sequence
	if ($fastaHash{$seq_name}) {
		my $seq = $fastaHash{$seq_name};
		printf ">$seq_name\t" . "-hmm_score $score\t" . "-no._domains $num_domains\t" . "\n";
		printf "$seq\n";
	}
#  	my $seq = $fastaHash{$seq_name};
#  	unless ("$seq" eq /^\s$/) { printf ">$seq_name\n$seq\n" };
}

exit;

