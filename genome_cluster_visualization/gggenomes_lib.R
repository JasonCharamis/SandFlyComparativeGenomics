# optionally install ggtree to plot genomes next to trees
# https://bioconductor.org/packages/release/bioc/html/ggtree.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

##Library to visualize genome clusters using gggenomes
package_install <- function ( package_name ) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  }
  
  else {
    print ( (sprintf("%s %s",package_name, "is not installed. Installing it!")))

    if ( package_name %in% BiocManager::available() ) {
      BiocManager::install(package_name)
    }
    
    else {
      install.packages(package_name)
    }
    
  }
}

dependencies <- c("tidyverse","dplyr","tibble","gggenomes")

for (i in dependencies) {
  package_install(i)
}



visualize_clusters <- function ( gff3_file, synteny=NULL, cluster_coord = NULL, ... ) {
  
  input_check_and_read <- function (input_file) {
    
    file_extension <- tools::file_ext(input_file)
    
    if ( grepl("bed", file_extension, ignore.case = TRUE) ) {
      gene_coordinates <- read_bed(input_file) ## if bed format is correct (strand is the last column, then gene orientation is automatically parsed)
    } 
    
    else if (grepl("gff3", file_extension, ignore.case = TRUE)) {
      gene_coordinates <- read_gff3(input_file)
      #gene_coordinates <- gene_coordinates %>% filter(grepl("gene|mRNA|exon|CDS",type))
    } 
    
    else { 
      stop("Input file format is not correct: Try bed or gff3 files") 
    }
    
  }
  
  ## reverse order chromosomes/species, default is to keep initial order as seen in the gff3 file
  genes <- input_check_and_read (gff3_file)
  
  ## get length of each chromosome instance and sort based on length, the largest goes down
  s0 <- genes %>% group_by(seq_id) %>% summarize(length=max(end)-min(start)) %>% arrange(length) 
  
  x <- full_join(genes, s0, by ='seq_id') ## include sequence length in tibble df
  x <- x %>% group_by(seq_id) %>% arrange(length) ## arrange genes based on sequence length
  
  if (is.null(synteny)) {
    #print ("Please provide a paf (default miniprot alignment output) file.")
    print ("Visualizing single clusters.")
    p <- gggenomes(genes = x, links = l0) 
    plot <- p + geom_seq() + geom_bin_label() + geom_gene()
  }
  
  else {
    
    links <- read_paf(synteny)  ## get synteny information from miniprot alignment
    l0 <- tibble(seq_id = links$seq_id2,
                 start = links$start2,
                 end = links$end2,
                 seq_id2 = links$seq_id,
                 start2 = links$start,
                 end2 = links$end,
                 hom = 1-links$de)
    
    p <- gggenomes(seqs = s0, genes = x, links = l0) 
    
    plot <- p +
      geom_seq() +         # draw contig/chromosome lines
      geom_bin_label() +  # draw contig/chromosome lines
      geom_gene() +        # draw genes as arrow
      geom_link(aes(fill=hom), color = "snow2", offset = 0.05) +  # draw some connections between syntenic regions
      scale_fill_gradient(name = "Sequence Homology",low = "snow1", high = "snow2") # color is based on sequence homology 
    
  }
  
  if ( exists("plot") == TRUE) {
    return (plot)
  }
}
