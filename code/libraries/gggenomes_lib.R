
#' gggenomes_lib
#'
#' @description
#' Library of functions to visualize genomic regions using the gggenomes R package (https://github.com/thackl/gggenomes).
#'
#' @import gggenomes
#' @importFrom gggenomes read_bed
#' @importFrom gggenomes read_gff3
#' @importFrom gggenomes geom_link
#' @importFrom gggenomes scale_fill_gradient
#' @importFrom dplyr filter
#' @importFrom tidyverse %>%
#' @importFrom tidyverse group_by
#' @importFrom tidyverse summarize
#' @importFrom tidyverse arrange
#' @importFrom tidyr full_join
#' @importFrom tools file_ext


#' .load_packages
#' @param tools A character vector specifying the packages to load.
#'
#' @return Nothing is returned; the function is used for loading packages.
#' 
#' @importFrom BiocManager install
#' @importFrom BiocManager available
#'
#' @export

.load_packages <- function( tools ) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version, ]
  
  for ( pkg in tools ) {
    if ( pkg %in% tmp$Package ) {
      library (pkg, character.only = TRUE)
    } else {
      print(sprintf("%s %s", pkg, "is not installed. Installing it!"))
      
      if ( pkg %in% BiocManager::available(pkg) ) {
        BiocManager::install(pkg, dependencies = TRUE, update = TRUE)
      } else {
        install.packages(pkg, dependencies = TRUE, ask = FALSE)
      }
    }
  }
}

dependencies <- c("tidyverse",
                  "dplyr",
                  "tibble",
                  "gggenomes")

.load_packages(dependencies)


#' @title .load_coordinate_file
#'
#' This function loads a coordinate file in BED, GFF3, or GTF/GFF2 format.
#'
#' @param coordinate_file Path to the coordinate file.
#'
#' @return An object representing the loaded coordinate file.
#'
#' @details The function checks the type of the coordinate file based on its extension 
#' (BED, GFF3, or GTF/GFF2) and reads it accordingly using appropriate functions 
#' (\code{read_bed}, \code{read_gff3}). If the provided file is not of the expected 
#' format, it throws an error message.
#'
#' @examples
#' # Load a BED file
#' .load_coordinate_file("example.bed")
#' # Load a GFF3 file
#' .load_coordinate_file("example.gff3")
#'
#' @export

.load_coordinate_file <- function( coordinate_file ) {
  if (typeof(coordinate_file) == "character") {
    if (file.exists(coordinate_file)) {
      if (any(grepl(".bed", coordinate_file))) {
        coordinate_file_obj <- gggenomes::read_bed(coordinate_file)
      } else if ( any(grepl(".gff3|.gff|.gtf", coordinate_file))) {
        coordinate_file_obj <- gggenomes::read_gff3(coordinate_file)
      } else {
        stop("Please provide a valid coordinate_file (BED, GFF3 or GTF/GFF2) containing the coordinates of the regions to be analyzed.")
      }
      
      return(coordinate_file_obj)
    }
  }
}


#' visualize_clusters
#' Main genome visualization function of the library
#' If multiple seq_ids are identified but the synteny argument is NULL ...
#  .. it will print each of the seq_ids independently
#  If the synteny argument is not NULL ...
#  .. it will also draw color-coded synteny based on the alignment information from the paf file
#' @param coordinate_file Path to the coordinate file.
#' @param synteny Path to the synteny file. Default is NULL.
#' @param cluster_coord Coordinates of  genomic regions. Default is NULL.
#' 
#' @return A ggplot object representing the visualization of genomic regions
#' @export

visualize_clusters <- function ( coordinate_file, 
                                 synteny=NULL, 
                                 cluster_coord = NULL,
                                 save = TRUE,
                                 output = NULL
                               ) {
  
  # Load coordinate file
  genes <- .load_coordinate_file (coordinate_file)
  
  # Get length of each chromosome instance and sort based on length, the largest goes down
  s0 <- genes %>% 
        group_by(seq_id) %>% 
        summarize( length = max(end)-min(start) ) %>% 
        arrange(length) 
  
  x <- full_join(genes, s0, by ='seq_id') # Include sequence length in tibble df
  
  x <- x %>% # Arrange genes based on sequence length in order to print the larger after the smaller ones
       group_by(seq_id) %>% 
       arrange(length) 
  
  if (is.null(synteny)) {
    print ("Visualizing single clusters.")
    
    p <- gggenomes(genes = x, links = l0) 
    
    plot <- p + 
            geom_seq() + # Draw contig/chromosome lines
            geom_bin_label() + # Print chromosome label
            geom_gene() # Draw genes as arrow
  } else {
      links <- read_paf(synteny)  # Get synteny information from minimap2 alignment
      
      l0 <- tibble(seq_id = links$seq_id2,
                   start = links$start2,
                   end = links$end2,
                   seq_id2 = links$seq_id,
                   start2 = links$start,
                   end2 = links$end,
                   hom = 1-links$de)
      
      p <- gggenomes(seqs = s0, 
                     genes = x, 
                     links = l0
                     ) 
      
      plot <- p +
        geom_seq() +         # Draw contig/chromosome lines
        geom_bin_label() +  # Print chromosome label
        geom_gene() +        # Draw genes as arrow
        geom_link(aes(fill=hom), color = "snow2", offset = 0.05) +  # Draw some connections between syntenic regions
        scale_fill_gradient(name = "Sequence Homology",low = "snow1", high = "snow2") # Color represents sequence homology 
  }
  
  
  # Export plot with provided options
  if (exists("plot")) {
    if (save == TRUE) {
      if (is.null(output)) {
        if (typeof(tree) == "character") {
          ggsave(plot = plot, sprintf("Genome_localization.svg", dpi = 600) )
          message <- "Plot saved as 'Genome_localization.svg'"
          print (message)
        } 
      } else {
        ggsave(plot = plot, output, dpi = 600)
        print(paste("Tree plotted and saved as", output))
      }
    } else {
      print("Plot will not be saved! Use the options save = TRUE and output = <OUTPUT_NAME> for saving the output plot.")
    }
    
    return ( plot )	
  }
}

