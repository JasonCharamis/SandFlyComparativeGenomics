## Perl script to combine the modified edgeR output with orthology information from Orthofinder output ##  
## Makes use of the 'Orthogroups/Orthogroups.txt' file ##

#use strict;
use warnings;

open ( IN, $ARGV[0] );
open ( IN1, $ARGV[1] );

my %DE_genes = ();
my %file = ();
my %geneids = ();
my $header = ();
my %orthologs = ();

## load DE gene IDs ##
while ( my $line = <IN> ) {
    chomp $line;
    if ( $line =~ /Geneid/ ) { $header = $line; next; }
    
    my @f = split (/\t/,$line);
    $DE_genes{$f[0]}=$line;
}

close (IN);

my %seen = ();
my @lines = ();

while ( my $line = <IN1> ) {

    chomp $line;
    push ( @lines, $line);
    my @f = split (/\t/,$line);
}

close (IN1);

## extract lines with orthologs for DE genes ##
my @h = ();

foreach ( keys %DE_genes ) {
    @h = grep(/$_/, @lines);
}


## associate DE genes with orthologs (if exist) ##
my @n = ();

my %ogs = ();
my %oggs = ();
my $gene = ();
my $n = ();
my $gen = ();
my @genes = ();

foreach ( @h ) {
    chomp $_;
    my @f = split (/:/,$_);
    ( my $ortho ) = ( $_ =~ /^(OG\w+):/g );
    ( my @number_of_genes ) = ( $f[1] =~ /\S+/g);   
    $n = scalar(@number_of_genes);
    $gen = $f[1];
    $gen =~ s/ /,/g;
    $gen =~ s/^,| //g;    
    @vc = split (/,/,$gen);
    for my $xc ( @vc ) {
        if ( $xc =~ /(aculy\w+)\.\d/ ) {
            $xc =~ s/\.\d//g;
            $ogs{$xc}="$ortho\_$n\_genes";
            $oggs{$xc}="$gen";
        }
    }
}


my %orthology = ();

foreach ( keys ( %oggs ) ) {

    if ( $oggs{$_} =~ /tetur/ ) {
        $orthology{$_}="Turticae";
    }
    elsif ( $oggs{$_} =~ /Dpter/ ) {
        $orthology{$_}="Dpter";
    }
    elsif ( $oggs{$_} =~ /ISCP/ ) {
        $orthology{$_}="Iscap";
    }
    elsif ( $oggs{$_} =~ /Mocc/ ) {
        $orthology{$_}="Mocc";
    }

    else { $orthology{$_}="Aculy_specific"; }
}

## add extra column in headers and print ##
my @hs = split (/\t/,$header);
print "$hs[0]\t$hs[1]\t$hs[2]\tOrthology\tOG_number\tOG_in_group\t";

for my $left ( 3..scalar(@hs)-1) {
    if ( $left <= scalar(@hs)-2 ) {
        print "$hs[$left]\t";
    }

    elsif ( $left == scalar(@hs)-1 ) {
        print "$hs[$left]\n";
    }   
}


# print new output with orthology 
foreach ( keys %DE_genes ) {

    my @j = split (/\t/,$DE_genes{$_});
    print "$j[0]\t$j[1]\t$j[2]\t$orthology{$_}\t$ogs{$_}\t$oggs{$_}\t";

    for my $lefti ( 3..scalar(@j)-1) {
        if ( $lefti <= scalar(@j) - 2) {
            print "$j[$lefti]\t";
        }
        elsif ( $lefti == scalar(@j) - 1 ) {
            print "$j[$lefti]\n";
        }
    }
}
