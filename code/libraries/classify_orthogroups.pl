
#use strict;
use warnings;

# Print script usage
sub print_usage {
    print <<"END_USAGE";

pwd: Orthofinder/Orthogroups
Usage: perl classify_orthogroups.pl Orthogroups.GeneCounts.tsv Orthogroups_Unassigned.tsv Orthogroups.txt

Options:
  -h, --help    Display this help message

END_USAGE
}

# Process command line arguments
if (@ARGV == 0 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help') {
    print_usage();
    exit;
}

# Scalars, arrays and hashes for parsing original file and counting number of species per OG
my @f = ();
my @count = ();
my $number_of_cols = ();
my %species_counter = ();

# Hashes for assigning OG id based on their presence
my %phlebotominae = ();
my %phlebotomus = ();
my %lutzomyia = ();
my %schw = ();
my %nematocera = ();
my %brachycera = ();

# Assign OG categories
my %universal = ();
my %universal_single_copy = ();
my %species_specific = ();
my %phlebotomus_specific = ();
my %lutzomyia_specific = ();
my %phlebotominae_wide = ();
my %phlebotominae_patchy = ();
my %nematocera_patchy = ();
my %diptera_patchy = ();

# Scalars, arrays and arrays for parsing the file, assigning OG id with species and saving the contents into a hash
my @file = ();
my %file = ();
my @h = ();
my @j = ();
my @y = ();
my $hg = ();

# Scalars, arrays and hashes for calculating the number of genes per category (universal, phlebotominae_wide etc) per species and printing results
my %data = ();

open ( OUT0, ">orthology_distribution_for_ggplot.txt" );
open ( OUT1, ">orthology_distribution.tsv" );

#=============================================================# PARSE THROUGH THE Orthogroups_Unassigned.tsv Orthogroups.txt FILES AND GET UNASSIGNED OG IDs #==================================================#

open ( IN2, $ARGV[1] ); # Orthogroups_UnassignedGenes.tsv

my %species_unassigned = (); # Hash to save the number of unassigned genes per species

while (my $line2 = <IN2>) {
    chomp $line2;
    
    if ($. == 1) {
        my @tabs = split(/\t/, $line2);
	
        for my $n (1 .. scalar(@tabs) - 1) {
            $species_unassigned{$tabs[$n]} = 1;
        }
    }
}


my %species_unassigned_counts = ();

# Load species unassigned OGs
foreach my $key (keys %species_unassigned) {
    my $grep_output = `grep -v "Orthogroup" $ARGV[1] | grep $key | cut -f1`;
    my @lines = split("\\n", $grep_output);
    my $nested_count = 0;

    foreach my $line (@lines) {
        $species_unassigned{$key}{$line} = 1;
        $nested_count++;
	$species_unassigned_counts{$key}=$nested_count;
    }
}


#=============================================================# PARSE THROUGH THE Orthogroups.GeneCount.tsv FILE, COUNT OG REPRESENTATION AND GENERATE RELEVANT ORTHOGROUP CATEGORIES #======================#

open ( IN, $ARGV[0] ); # Orthogroups.GeneCount.tsv

while ( my $line = <IN> ) {
    
    chomp $line;
    
    if ( $line =~ /\w/ ) {
	push ( @file, $line);
    }

    # First characterize orthogroups based on their presence (universal, widespread, phlebotominae, culicidae, chironomidae, nematocera, brachycera, none - note that this is per species )    
    if ( $. > 1 ) { 

	@f = split (/\t/,$line); 
	$number_of_cols = scalar(@f);
    
	my $i = 0;

	for my $species (1..$number_of_cols-2) { 	
	    if ( $f[$species] > 0 ) { 
		$i++;
		$species_counter{$f[0]} = $i;
		
		if ( $f[$species] == 1 ) {
		    $file{$f[0]}{$species} = $f[$species];
		}
	    }
	}

	my $p = 0;
	for my $phlebotomines (9..$number_of_cols-2) { # OGs present in phlebotominae with maximum species counts
	    if ( $f[$phlebotomines] > 0 ) {
		$p++;
		$phlebotominae{$f[0]} = $p;
	    }
	}

	my $ph = 0;
	for my $phlebotomus (9..11, 14..16, 18, 19) {
	    if ( $f[$phlebotomus] > 0 ) {
		$ph++;
		$phlebotomus{$f[0]} = $ph;
	}

	my $ltz = 0;
	for my $lutzomyia (12, 13) {
	    if ( $f[$lutzomyia] > 0 ) {
		$ltz++;
		$lutzomyia{$f[0]} = $ltz;
	    }
	}

	my $stz = 0;
	for my $schw (17) {
	    if ( $f[$schw] > 0 ) {
		$stz++;
		$schw{$f[0]} = $stz;
	    }
	}	

	my $m = 0;
	for my $mosquitoes_midge (1,2,3,4,5,8) {
	    if ( $f[$mosquitoes_midge] > 0 ) {
		$m++;
		$nematocera{$f[0]} = $m;
	    }
	}

	my $d = 0;
	for my $flies (6,7) {
	    if ( $f[$flies] > 0 ) {
		$d++;
		$brachycera{$f[0]} = $d;
	}
      }
    }
  }
}


my $universal_single_copy = 0;
my $og = ();
    
for $og ( keys %file ) {
    $universal_single_copy = 0;
    
    foreach ( keys %{ $file{$og} } ) { # Universal single-copy
	$universal_single_copy += $file{$og}{$_} ;
    }

    if ( $universal_single_copy ==  ( $number_of_cols-2 ) ) {
	$universal_single_copy{$og} = 0;
    }
}
 
foreach ( sort keys %species_counter) {  # Universal: orthogroups found in all or all-but-two species
    if ( $species_counter{$_} >= $number_of_cols-4 ) {
	unless ( exists $universal_single_copy{$_} ) {
	    $universal{$_} = $species_counter{$_};
	}
    }
    
    elsif ( $species_counter{$_} < 2 ) {
	$species_specific{$_} = $species_counter{$_};
    }
}

foreach ( sort keys %phlebotominae ) { # Widespread: orthogroups found in all or all-but-two Phlebotominae species
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $brachycera{$_} ) {
		    unless ( exists $nematocera{$_} ) {
			if ( $phlebotominae{$_} >= 9 ) {
			    $phlebotominae_wide{$_} = $phlebotominae{$_};
			}
		    }
		}
	    }
	}
    }
}


foreach ( sort keys %phlebotomus ) { # Phlebotomus-specific: orthogroups found only in Phlebotomus and in all or all-but-one species
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $brachycera{$_} ) {
		    unless ( exists $nematocera{$_} ) {
			unless ( exists $phlebotominae_wide{$_} ) {
			    unless ( exists $lutzomyia{$_} ) {
				unless ( exists $schw{$_} ) {			       
				    if ( exists $phlebotomus{$_} && $phlebotomus{$_} >= 7 ) {
					$phlebotomus_specific{$_} = $phlebotomus{$_};
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

foreach ( sort keys %lutzomyia ) { # Lutzomyia-specific
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $brachycera{$_} ) {
		    unless ( exists $nematocera{$_} ) {
			unless ( exists $phlebotominae_wide{$_} ) {
			    unless ( exists $phlebotominae_patchy{$_} ) {
				unless ( exists $phlebotomus{$_} ) {
				    unless ( exists $schw{$_} ) {			       
					if ( exists $lutzomyia{$_} && $lutzomyia{$_} == 2 ) {
					    $lutzomyia_specific{$_} = $lutzomyia{$_};
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

foreach ( sort keys %phlebotominae ) { # Phlebotominae-patchy: orthogroups present only in phlebotomines with patchy distribution (general category)
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $brachycera{$_} ) {
		    unless ( exists $nematocera{$_} ) {
			unless ( exists $phlebotominae_wide{$_} ) {
			    unless ( exists $phlebotomus_specific{$_} ) {
				unless ( exists $lutzomyia_specific{$_} ) {
				    if ( $phlebotominae{$_} >= 2 && $phlebotominae{$_} <= 8 ) {
					$phlebotominae_patchy{$_} = $phlebotominae{$_};
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

foreach ( sort keys %nematocera) { # Nematocera
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $brachycera{$_} ) {
		    unless ( exists $phlebotominae_wide{$_} ) {
			unless ( exists $phlebotominae_patchy{$_} ) {
			    unless ( exists $phlebotomus_specific{$_} ) {
				unless ( exists $lutzomyia_specific{$_} ) {
				    if (  $species_counter{$_} >= 2  && $species_counter{$_} < 19) {
					$nematocera_patchy{$_} = $nematocera{$_};
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


foreach ( sort keys %brachycera) { # Diptera
    unless ( exists $universal_single_copy{$_} ) {
	unless ( exists $universal{$_} ) {
	    unless ( exists $species_specific{$_} ) {
		unless ( exists $nematocera_patchy{$_} ) {
		    unless ( exists $phlebotominae_wide{$_} ) {
			unless ( exists $phlebotominae_patchy{$_} ) {
			    unless ( exists $phlebotomus_specific{$_} ) {
				unless ( exists $lutzomyia_specific{$_} ) {
				    $diptera_patchy{$_} = $brachycera{$_};			      
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


#=============================================================# PARSE THROUGH THE FILE THE OBTAIN THE OGS FOR EACH CATEGORY PER SPECIES AS WELL AS THE GENE NUMBER OF EACH SUCH INSTANCE #===================#

for my $nl (0..scalar(@file)-1) {

    # Create two arrays: the first one will parse the file line-by-line, while the second one tab-by-tab
    my @j = split (/\n/,$file[0]);
    my @y = split (/\t/,$j[0] ); 

    if ($nl == 0) {
        next;
    }

    my @h = split (/\t/,$file[$nl]);

    my $orthogroup = $h[0];

    my $i = 0;
 
    for my $col (1..scalar(@h)-2) {
        my $species = $y[$col];
        my $gene_count = $h[$col];
        $data{$species}{$orthogroup} = $gene_count;
    }
}


my %gene_counts;

# Push OGs to specified categories
for my $sp (sort keys %data) {
    for my $og (sort keys %{$data{$sp}}) { 
        if (exists($universal_single_copy{$og})) {
            push(@{$data{$sp}{'Universal_single_copy'}}, $og);
        } elsif (exists($universal{$og})) {
            push(@{$data{$sp}{'Universal'}}, $og);
        } elsif (exists($phlebotominae_wide{$og})) {
            push(@{$data{$sp}{'Phlebotominae_wide'}}, $og);
        } elsif (exists($phlebotomus_specific{$og})) {
            push(@{$data{$sp}{'Phlebotomus_specific'}}, $og);
        } elsif (exists($lutzomyia_specific{$og})) {
            push(@{$data{$sp}{'Lutzomyia_specific'}}, $og);
        } elsif (exists($phlebotominae_patchy{$og})) {
            push(@{$data{$sp}{'Phlebotominae_patchy'}}, $og);
        } elsif (exists($nematocera_patchy{$og})) {
            push(@{$data{$sp}{'Nematocera_patchy'}}, $og);
        } elsif (exists($diptera_patchy{$og})) {
            push(@{$data{$sp}{'Diptera_patchy'}}, $og);
        } elsif (exists($species_specific{$og})) {
            push(@{$data{$sp}{'Species_specific'}}, $og);
	}
    }
}

my %seen = ();

for my $species (sort keys %data) {
    for my $array (sort keys %{$data{$species}}) {
        $gene_counts{$species}{$array} = 0; # Initialize the count to 0 for each category

        for my $og (@{$data{$species}{$array}}) {
            $gene_counts{$species}{$array} += $data{$species}{$og};
        }
    }
}


# Add the number of species-unassigned genes in the Species-specific subset
foreach my $species (keys %species_unassigned) {
    $gene_counts{$species}{'Species_specific'} += $species_unassigned_counts{$species};
}

my @categories = qw(Species Universal_single_copy Universal Phlebotominae_wide Phlebotomus_specific Lutzomyia_specific Phlebotominae_patchy Nematocera_patchy Diptera_patchy Species_specific); # Define custom order

# Generate headers for the two files with common structure ( Species as Col1 and Relevant Categories in the rest of the columns )
my $printed = 0; # Flag to keep track of whether OUT0 has been printed

if (!$printed) { # Print header ONCE, with custom order

    # Print the header
    print OUT0 "Species\tNumber_of_Genes\tType\n";

    my $num_keys = scalar ( @categories );
    my $count = 0;

    for my $orthology (@categories) {
	$count++;
		
	if ($count < $num_keys) {
	    print OUT1 "$orthology\t";
	} else {
	    print OUT1 "$orthology\n";
	}
    }
    
    $printed = 1; # Set the flag to indicate that printing has been done
}


# Print gene counts per OG per species
for my $species ( sort keys %gene_counts) {
    
    print OUT1 "$species\t"; # First print the sorted species names

    my $num_orthology = scalar ( @categories );
    my $orthology_count = 0;
    
    for my $orthology (@categories) {
	$orthology_count++;

	my %seen;
	
	unless (exists $seen{$species}) {
	    unless ( $orthology eq 'Species' ) {
		print OUT0 "$species\t$gene_counts{$species}{$orthology}\t$orthology\n";
	    }
	    
	    print OUT1 "$gene_counts{$species}{$orthology}";
	    $seen{$species}{$orthology} = 1;
	}

	if ($orthology_count < $num_orthology) {
	    print OUT1 "\t";
	} else {
	    print OUT1 "\n";
	}
    }
}
