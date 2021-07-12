#!/usr/bin/perl
use strict;
my %geneCover;
my %geneLength;
open IN, "$ARGV[0]";
while(<IN>){
	my @l = split /\s+/, $_;
	if($l[3] >= 2){
		$geneCover{$l[7]}{$l[1]} = $l[2];
	}else{
		$geneCover{$l[7]}{0} = 0;
	}
	$geneLength{$l[7]} = $l[6] - $l[5] + 1;
}
close IN;

print"Gene\tGeneLength\tGeneCover\tGeneCoverRate\tPAV\n";
foreach my $gene (sort {$a cmp $b} keys %geneCover){
	my $depthLen = 0;
	foreach my $pos (sort {$a <=> $b} keys %{$geneCover{$gene}}){
	$depthLen = $depthLen + $geneCover{$gene}{$pos} - $pos;
	}
	$depthLen = $depthLen + 1;
	my $cover = $geneLength{$gene} > 0 ? $depthLen/$geneLength{$gene} : 0;
	my $pav;
	if($cover>=0.4){
		$pav = 1;
	}else{
		$pav = 0;
	}
	print"$gene\t$geneLength{$gene}\t$depthLen\t$cover\t$pav\n"; 
}
