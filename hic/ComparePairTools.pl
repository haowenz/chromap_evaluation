#!/usr/bin/env perl

use strict;
use warnings;

die "usage: a.pl a.pairs b.pairs\n" if (@ARGV == 0) ;

my %coordinate;
my %readUsed ;
open FP, $ARGV[0];
while (<FP>)
{
	next if (/^#/) ;
	chomp ;
	my @cols = split ;
	next if ($cols[1] eq "!" || $cols[3] eq "!") ;
	@{$coordinate{$cols[0]}} = ($cols[1], $cols[2], $cols[3], $cols[4]) ;
	$readUsed{$cols[0]} = 0;
}
close FP ;

my $aonly = 0;
my $bonly = 0;
my $intersect = 0;
my $intersectOneside = 0;
my $different = 0;
open FP, $ARGV[1] ;
while (<FP>)
{
	next if (/^#/) ;
	chomp ;
	my @cols = split ;
	next if ($cols[1] eq "!" || $cols[3] eq "!") ;
	if (!defined $coordinate{$cols[0]}) 
	{
		#print($cols[0], "\n");
		++$bonly ;
		next;
	}
	++$readUsed{$cols[0]};
	my $match = 0;
	my $width = 20;
	my $i;
	for ($i = 1 ; $i <= 3 ; $i += 2)
	{
		my $k = 4 - $i;
		my @ac = @{$coordinate{$cols[0]}};
		if ($cols[$i] ne $ac[0] && $cols[$i] ne $ac[2])
		{
			next ;
		}
		if ($cols[$i] eq $ac[0] && abs($cols[$i + 1] - $ac[1]) <= $width)
		{
			$match |= int(($i + 1)/2) ;
		}
		if ($cols[$i] eq $ac[2] && abs($cols[$i + 1] - $ac[3]) <= $width)
		{
			$match |= int(($i + 1)/2) ;
		}
	}
	
	if ($match == 0)
	{
		++$different ;
	}
	elsif ($match == 3)
	{
		++$intersect ;
	}
	else
	{
		#print($cols[0], "\n");
		++$intersectOneside; 
	}
}

for my $key (keys %readUsed) 
{
	if ($readUsed{$key} == 0)	
	{
		#print($key, "\n");
		++$aonly ;
	}
}

print(join("\t", ($aonly, $bonly, $different, $intersectOneside, $intersect)), "\n");
