use strict;
use warnings;

my %genes = ();
my @lines = ();

open ( IN, $ARGV[0] );

while ( my $line = <IN> ) {
    chomp $line;

    $line =~ s/PSKW_//g;

    #print header based on edgeR file architecture
    if ( $line =~ /Geneid/ ) {

        my @f = split (/\t/,$line);

        ## keep only geneid, log2FC, p-value and FDR ##
        print "Geneid\tFold_Change\tAnnotation\tP_value\tFDR\t";

        ## print raw read counts per replicate ##
        for my $p ( 7..scalar(@f) - 1 ) {

            ## no samples are removed by default ##
            #           unless ( $p =~ /\b9|13|14|17/ ) {
            if ( $p <= scalar(@f)-2 ) {
                print "$f[$p]\t";
            }
#       }
            elsif ( $p == scalar(@f-1)) {
                print "$f[$p]\n";
            }
        }
    }

    push ( @lines, $line);

}

open ( IN1, $ARGV[1] );

while ( my $line = <IN1> ) {
    chomp $line;
    $line =~ s/\.\d//g;
    my @f = split (/\t/,$line);

    ## make a hash with geneid - description associations ##
    unless ( exists $genes{$f[0]} ) {
        $genes{$f[0]} = $f[6];
    }

}


foreach ( @lines ) {
  ## combine edgeR output and annotation file to print the final output with
  ## only the desired columns from edgeR output and the gene description/annotation
   my @h = split (/\t/,$_);
   if ( exists $genes{ $h[0] } ) {
       print "$h[0]\t",dec(fold($h[3])),"\t","$genes{$h[0]}\t",scient($h[5]),"\t",scient($h[6]),"\t";

        for my $i ( 7..scalar(@h) - 1 ) {
        ## no samples are removed by default ##
            #unless ( $i =~ /\b9|13|14|17/ ) {

            if ( $i <= scalar(@h)-2 ) {
                print "$h[$i]\t";
            }   
#       }                                                                                                                                                                                                         
            elsif ( $i == scalar(@h-1)) {
                print "$h[$i]\n";
            }
        }
    }
}

sub dec {
    my @input = @_;
    my $line = $input[0];
    my $rounded = sprintf("%.3f", $line);
    return $rounded;   
}

sub scient {
    my @input = @_;
    my $line = $input[0];
    my $rounded = sprintf("%.2e", $line);
    return $rounded;   
}

sub fold {
    my @input = @_;
    my $line = $input[0];
    my $fold = ();
    
    if ( $line > 0 ) {
        $fold = 2**($line);
    }

    if ( $line < 0 ) {
        $fold = -2**abs(($line));
    }
    return $fold;
}
