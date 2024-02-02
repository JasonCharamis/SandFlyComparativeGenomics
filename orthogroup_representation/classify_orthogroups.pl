
use strict;
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

#scalars, arrays and hashes for parsing original file and counting number of species per OG
my @f = ();
my @count = ();
my $number_of_cols = ();
my %species_counter = ();

#hashes for assigning OG id based on their presence
my %phlebotominae = ();
my %phlebotomus = ();
my %lutzomyia = ();
my %schw = ();
my %nematocera = ();
my %brachycera = ();

#assign OG categories
my %universal = ();
my %universal_single_copy = ();
my %species_specific = ();
my %phlebotomus_specific = ();
my %lutzomyia_specific = ();
my %phlebotominae_wide = ();
my %phlebotominae_patchy = ();
my %nematocera_patchy = ();
my %diptera_patchy = ();

#scalars, arrays and arrays for parsing the file, assigning OG id with species and saving the contents into a hash
my @file = ();
my %file = ();
my @h = ();
my @j = ();
my @y = ();
my $hg = ();

#scalars, arrays and hashes for calculating the number of genes per category (universal, phlebotominae_wide etc) per species and printing results
my %data = ();

open ( IN, $ARGV[0] ); ## Orthogroups.GeneCount.tsv
open ( IN2, $ARGV[1] ); ## Orthogroups_UnassignedGenes.tsv
open ( IN3, $ARGV[2] ); ## Orthogroups.txt

open ( OUT0, ">orthology_distribution_for_R.txt" );
open ( OUT1, ">OGs_per_category_per_species.txt" );
open ( OUT2, ">orthology_distribution.tsv" );

##=============================================================## PARSE THROUGH THE Orthogroups_Unassigned.tsv Orthogroups.txt FILES AND GET UNASSIGNED OG IDs ##==================================================##

my %species_unassigned = (); ## hash to save the number of unassigned genes per species

while (my $line2 = <IN2>) {  
    if ($. == 1) {
        my @tabs = split(/\t/, $line2);
        for my $n (1 .. scalar(@tabs) - 2) {
            $species_unassigned{$tabs[$n]} = 1;
        }
    }
}

foreach (keys %species_unassigned) {
    $species_unassigned{$_} = `grep $_ $ARGV[1] | grep -v "Orthogroup" | cut -f1`;
}

my %OG2genes = (); ## hash to save the number of unassigned genes per species

	while (my $line3 = <IN3>) {
		chomp $line3;
		my @t = split(/: /, $line3);
		
		foreach ( @t ) {
		$OG2genes{$t[0]} = $t[1]; # OG to genes
		}
	}

##=============================================================## PARSE THROUGH THE Orthogroups.GeneCount.tsv FILE, COUNT OG REPRESENTATION AND GENERATE RELEVANT ORTHOGROUP CATEGORIES ##======================##
while ( my $line = <IN> ) {
    
    chomp $line;
    
    if ( $line =~ /\w/ ) {
	push ( @file, $line);
    }
  
    #first characterize orthogroups based on their presence (universal, widespread, phlebotominae, culicidae, chironomidae, nematocera, brachycera, none - note that this is per species )    
    if ( $. > 1 ) { ## skip header line with species names

	@f = split (/\t/,$line); #split file by columns
	$number_of_cols = scalar(@f);
    
	my $i = 0; #create a variable to use as counter when finding one of the species

	for my $species (1..$number_of_cols-2) { 	#skip orthogroup name and total (0 and 20)
	    if ( $f[$species] > 0 ) { #count number of species with higher than zero number of genes per orthogroup
		$i++;
		$species_counter{$f[0]} = $i;
		
		if ( $f[$species] == 1 ) {
		    $file{$f[0]}{$species} = $f[$species];
		}
    }
}

	my $p = 0;
	for my $phlebotomines (9..$number_of_cols-2) { #OGs present in phlebotominae with maximum species counts
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


my $universal_single_copy = 0;
my $og = ();
    
for $og ( keys %file ) {
    $universal_single_copy = 0;
    
    foreach ( keys %{ $file{$og} } ) {
	$universal_single_copy += $file{$og}{$_} ;
    }

    if ( $universal_single_copy ==  ( $number_of_cols-2 ) ) {
	$universal_single_copy{$og} = 0;
    }
}
 
foreach ( sort keys %species_counter) {  #if orthogroup is found in all or all-but-two species save the ids into the universal dictionary    
    if ( $species_counter{$_} >= $number_of_cols-4 ) {
	unless ( exists $universal_single_copy{$_} ) {
	    $universal{$_} = $species_counter{$_};
	}
    }
    
    elsif ( $species_counter{$_} < 2 ) {
	$species_specific{$_} = $species_counter{$_};
    }
}

foreach ( sort keys %phlebotominae ) { #orthogroups widespread (all or all-but-two) in phlebotominae    
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


foreach ( sort keys %phlebotomus ) { #orthogroups specific (all or all-but-one) to phlebotomus genus 
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

foreach ( sort keys %lutzomyia ) { #orthogroups specific (all) to lutzomyia genus 
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

foreach ( sort keys %phlebotominae ) { #orthogroups present only in phlebotomines, and with patchy distribution
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

foreach ( sort keys %nematocera) { #orthogroups present in nematocera, with patchy distribution
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


foreach ( sort keys %brachycera) { #orthogroups present in diptera, with patchy distribution
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



##=============================================================## PARSE THROUGH THE FILE THE OBTAIN THE OGS FOR EACH CATEGORY PER SPECIES AS WELL AS THE GENE NUMBER OF EACH SUCH INSTANCE ##===================##

for my $nl (0..scalar(@file)-1) {
    @j = split (/\n/,$file[0]);
    @y = split (/\t/,$j[0] ); # create an array with the species names

    if ($nl == 0) {
        next;
    }

    @h = split (/\t/,$file[$nl]);

    my $orthogroup = $h[0];

    for my $col (1..scalar(@h)-2) {
        my $species = $y[$col];
        my $gene_count = $h[$col];

        $data{$species}{$orthogroup} = $gene_count;
    }
}

my %subsets;
my %gene_counts;

for my $sp (sort keys %data) {
    for my $og (sort keys %{$data{$sp}}) { 
        if (exists($universal_single_copy{$og})) {
            push(@{$subsets{$sp}{Universal_single_copy}}, $og);
        } elsif (exists($universal{$og})) {
            push(@{$subsets{$sp}{Universal}}, $og);
        } elsif (exists($phlebotominae_wide{$og})) {
            push(@{$subsets{$sp}{Phlebotominae_wide}}, $og);
        } elsif (exists($phlebotomus_specific{$og})) {
            push(@{$subsets{$sp}{Phlebotomus_specific}}, $og);
        } elsif (exists($lutzomyia_specific{$og})) {
            push(@{$subsets{$sp}{Lutzomyia_specific}}, $og);
        } elsif (exists($phlebotominae_patchy{$og})) {
            push(@{$subsets{$sp}{Phlebotominae_patchy}}, $og);
        } elsif (exists($nematocera_patchy{$og})) {
            push(@{$subsets{$sp}{Nematocera_patchy}}, $og);
        } elsif (exists($diptera_patchy{$og})) {
            push(@{$subsets{$sp}{Diptera_patchy}}, $og);
        } elsif (exists($species_specific{$og})) {
            push(@{$subsets{$sp}{Species_specific}}, $og);
	}
    }

    my @unassigned_ogs = split (/\n/,$species_unassigned{$sp});
    foreach ( @unassigned_ogs ) {
	push(@{$subsets{$sp}{Species_specific}}, $_); ## Add the unassigned genes per species, as species-specific
    }
}


# Generate counts for each category per species
for my $sp (sort keys %subsets) {
    for my $array (sort keys %{$subsets{$sp}}) {
	for my $og (@{$subsets{$sp}{$array}}) {
            if (exists $data{$sp}{$og}) {
                $gene_counts{$sp}{$array} += $data{$sp}{$og};	
            }
        }
    }
}

print OUT0 "Species\tNumber_of_Genes\tType\n";

foreach my $sp ( sort keys %gene_counts ) {
    for my $category ( keys %{$gene_counts{$sp}}) {
        print OUT0 "$sp\t$gene_counts{$sp}{$category}\t$category\n";
    }
}


my @categories = qw(Species Universal_single_copy Universal Phlebotominae_wide Phlebotomus_specific Lutzomyia_specific Phlebotominae_patchy Nematocera_patchy Diptera_patchy Species_specific); ## Define custom order

## Generate headers for the two files with common structure ( Species as Col1 and Relevant Categories in the rest of the columns )
my $printed = 0; # Flag to keep track of whether OUT0 has been printed

if (!$printed) { ## Print header ONCE, with custom order

    my $num_keys = scalar ( @categories );
    my $count = 0;

    for my $orthology (@categories) {
	$count++;
		
	if ($count < $num_keys) {
	    print OUT1 "$orthology\t";
	    print OUT2 "$orthology\t";
	} else {
	    print OUT1 "$orthology\n";
	    print OUT2 "$orthology\n";
	}
    }
    
    $printed = 1; # Set the flag to indicate that printing has been done
}


## Print (OUT1) OG IDs per Category per species and (OUT2) actual output of gene counting per OG per species
for my $sp ( sort keys %subsets) { 
    print OUT1 "$sp\t"; ## First print the sorted species names
    print OUT2 "$sp\t"; ## First print the sorted species names

    my $num_orthology = scalar ( @categories );
    my $orthology_count = 0;
    
    for my $orthology (@categories) {
	$orthology_count++;

	my %seen;
	my %seen2;

	if (defined($subsets{$sp}{$orthology}) && ref($subsets{$sp}{$orthology}) eq 'ARRAY') {
	   # print OUT1 join(",", @{$subsets{$sp}{$orthology}});

	    for my $og (@{$subsets{$sp}{$orthology}}) {
		if (exists $OG2genes{$og}) {
		    unless (exists $seen{$og}) {
			print OUT1 "$OG2genes{$og}" . ",";
			$seen{$og} = 1;
		    }
		}
	    }
	} else {
	    print("Cannot read this line.\n");
	}

	unless (exists $seen2{$sp}) {
	    print OUT2 "$gene_counts{$sp}{$orthology}";
	    $seen2{$sp}{$orthology} = 1;
	}

	if ($orthology_count < $num_orthology) {
	    print OUT1 "\t";
	    print OUT2 "\t";
	} else {
	    print OUT1 "\n";
	    print OUT2 "\n";
	}
    }
}
