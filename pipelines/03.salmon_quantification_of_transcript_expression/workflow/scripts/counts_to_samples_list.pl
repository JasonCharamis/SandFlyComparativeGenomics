
use strict;
use warnings;

open ( IN, $ARGV[0] );

while ( my $line = <IN> ) {
    unless ( $line =~ /\w/ ) {
	next;
    }

    if ( $. == 1 ) {   
	chomp $line;
    	my @f = split (/\t/,$line);

	for my $i (1..scalar @f) { 
	    my $sample = $f[$i] ; 
	    $sample =~ s/\d+$|\_\d+$//g; 
	    print "$f[$i]\t$sample\n";
	}
    }
}
