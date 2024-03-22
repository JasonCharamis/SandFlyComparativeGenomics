
#' @title .load_packages
#' @description This function checks if a package is installed. If not, it installs the package using BiocManager if available, otherwise using install.packages.
#' @param tools A character vector of package names to be checked and installed.
#' @return NULL
#' @export

# Function to check if a package is installed, and if not, install it.
# Function to install or load a package

.load_packages <- function(tools) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built,
                                      1,
                                      1)))
  tmp <- tmp[as.numeric(substr(tmp$Built,
                                1,
                                1)) == max_version,
             ]

  for (pkg in tools) {
    if (pkg %in% tmp$Package) {
      suppressPackageStartupMessages(library(pkg,
                                             character.only = TRUE))
    } else {
      print(sprintf("%s %s",
                    pkg,
                    "is not installed. Installing it!"))
      
      if (pkg %in% BiocManager::available(pkg)) {
        BiocManager::install(pkg,
                             dependencies = TRUE,
                             update = TRUE)
      } else {
        install.packages(pkg,
                         dependencies = TRUE,
                         ask = FALSE)
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c( "ggplot2",
                   "dplyr",
                   "scales",
                   "tidyverse")

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

#' @title count_genes_sliding_window
#'
#' This function counts genes from a BED, GFF3, or GTF/GFF2 file using a sliding window approach and generates a plot showing gene counts across the genome.
#'
#' @param input_file Path to the input BED file.
#' @param sliding_window_size The size of the sliding window for counting genes (default is 50000).
#' @param bar_color Color for the bars in the plot (default is "black").
#' @param title Title for the plot (default is constructed based on the sliding window size).
#' @param save Option to save plot (default is FALSE).
#'
#' @return A ggplot object representing gene counts across the genome.
#'
#' @details This function first loads the coordinate file using the \code{.load_coordinate_file} function, then plots the gene counts per sliding window across the chromosomes provided in the input BED file. The plot is generated using ggplot2 package.
#'
#' @examples
#' # Count genes from a BED file with default parameters
#' count_genes_from_bed("example.bed")
#' # Count genes from a BED file with custom sliding window size and title
#' count_genes_from_bed("example.bed", sliding_window_size = 10000, title = "Gene Counts with 10kb Window Size")
#'
#' @export


count_genes_sliding_window <- function( coordinate_file,
                                        sliding_window_size = 50000,
                                        bar_color = "black",
                                        title = paste("Gene Counts using a", sliding_window_size, "window size", collapse = " "),
                                        save = FALSE ) {

			
		     genes <- .load_coordinate_file( coordinate_file )

		     # Plot the gene counts per sliding window across the provided chromosomes
		     
		     plot <- ggplot( genes, aes( x = start, y = after_stat(count)), fill = bar_color) +		     
		     	       geom_histogram( binwidth = sliding_window_size, color = "black", alpha = 0.7, position = "identity" ) +
      			     labs( title = title, x = "Genome Position", y = "Gene Count" ) +
      			     facet_wrap(~ seq_id, scales = "free_x", ncol = length(unique(genes$seq_id))) +
      			     theme_minimal() +
      			     theme( plot.title = element_text(hjust = 0.5, face = "bold"),
      				   axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
      				   axis.text.y = element_text(face = "bold") )
      					  

		     # Export plot with provided options
		     
		     if (exists("plot")) {
		     	if (save == TRUE) {
		     	   if (is.null(output)) {
		     	      if (typeof(tree) == "character") {
		     	      	 ggsave(plot = plot, sprintf("Gene_counts_across_genome.svg", dpi = 600) )
		     		       message <- "Plot saved as 'Gene_counts_across_genome.svg'"
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

