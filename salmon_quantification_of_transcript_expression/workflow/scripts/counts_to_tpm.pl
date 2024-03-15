
# Author: Panagiotis Ioannidis

#!/usr/bin/perl -w

use strict;
use warnings;

if( scalar ( @ARGV ) != 1 ){
	my @temp=split(/\//, $0);
	die
"

USAGE: $temp[-1] <featureCounts output>


";
}
# This script will take in the featureCounts output 
# and calculate TPM expression values for each gene

my $cf = 6; # counts' field: fields from this number onward contain read count values (zero-based)

my $lf = 5; # length field; the field number containing the gene length (cumulative exon length)
 
my %output = ();

my @sum_rpk = ();

open ( IN, $ARGV[0] ) or die;

while ( my $line = <IN> ) {
	unless ( $line =~ /^\S/ ) { next; }
	
	if ( $line =~ /^#/ ) {
		#print $line;
		next;
	}

	chomp ($line);
	my @f = split (/\t/, $line);
	
	if ( $line =~ /^Geneid/ ) { # if the line is the header row, then only print the fields you're interested in
		print "$f[0]"; # print the header of the first field (gene name)
		for ( my $i = $cf; $i < scalar(@f); $i++ ) {
			$f[$i] =~ s/\/out\.bam//;
			print "\t$f[$i]";
		}
		print "\n";
		next;
	}
	
	my $gene_len = $f[$lf] / 1000; #length of gene in kbp

	my @tmp = (); # will hold the RPK values for each sample/replicate of a gene
	
	for ( my $i = $cf; $i < scalar(@f); $i++ ) {
		$f[$i] /= $gene_len;
		push ( @tmp, $f[$i] );
		
		$sum_rpk[ $i - $cf ] += $f[$i]; # increase the sum of RPKs for each sample/replicate
	}

	$output{ $f[0] } = [ @tmp ]; # assign the list of RPKs to the gene

	#print STDERR "$f[0]";
	#my @tmp2 = @{ $output{ $f[0] } };
	#foreach my $tmp2 (@tmp2) {
	#	print STDERR "\t$tmp2";
	#}
	#print STDERR "\n";
}

close (IN);

for ( my $i = 0; $i < scalar( @sum_rpk ); $i++ ) {
	$sum_rpk[$i] /= 1000000; # divide each sum by 1,000,000
}

#print STDERR join ( "\t", @sum_rpk ), "\n";

my @output = keys ( %output );

#foreach my $tmp (@output) {
#	my ( $tmp1 ) = ( $tmp =~ /g(\d+)/ );
#	print STDERR "$tmp\t$tmp1\n";
#}

#@output = sort {
#	my ( $aVal ) = ( $a =~ /LOC(\d+)/ );
#	my ( $bVal ) = ( $b =~ /LOC(\d+)/ );
#	$aVal <=> $bVal;
#} @output;

@output = sort (@output);

foreach my $gene ( @output ) {
	print $gene;
	my @rpk = @{ $output{ $gene } }; # de-reference the array containing the RPK values
	
	#print STDERR $gene;
	#foreach my $tmp (@rpk) {
	#	print STDERR "\t$tmp";
	#}
	#print STDERR "\n";
	
	for ( my $k = 0; $k < scalar(@rpk); $k++ ) {
		my $tpm = $rpk[$k] / $sum_rpk[$k];
		$tpm = sprintf ( "%.3f", $tpm );
		
		print "\t$tpm";
	}
	print "\n";
}

