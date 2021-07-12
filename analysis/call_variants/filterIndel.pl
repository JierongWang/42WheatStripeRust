#!/usr/bin/perl
use strict;
open IN, "$ARGV[0]";
while(<IN>){
	chomp;
	if(/^#CHROM/){
		print"$_\n";
	}
	next if(/^#/);
	my @l = split /\s+/, $_;
	next if($l[4]=~/,/);
	my $refLen = length($l[3]);
	my $altLen = length($l[4]);
	if($refLen <= 50 && $altLen <= 50){
		print"$_\n";
	}
}
close IN;
