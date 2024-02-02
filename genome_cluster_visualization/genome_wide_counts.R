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

dependencies <- c("gggenes","tidyverse")

for (i in dependencies) {
  package_install(i)
}



### LOCALIZATION OF P450 GENES IN LUTZOMYIA LONGIPALPIS GENOME ###

genes <- read.delim("Phlebotomus_papatasi.P450s.bed2")
genes <- read.delim("Lutzomyia_longipalpis.P450s.manually_curated.bed")
#genes <- genes %>% mutate(seq_id = gsub(".*chr", "chr", seq_id)) %>% mutate (seq_id = gsub("_"," ",seq_id))
#genes %>% filter (type=="gene")

ggplot(genes, aes(xmin = start, xmax = end, y = seq_id)) +
  geom_gene_arrow() +
  theme_genes()

# Plot the gene counts
# Assuming your bin_width is defined
bin_width <- 50000
color <- "black"

p450_localization <- ggplot(genes, aes(x = start, y = after_stat(count), fill = color)) +
  geom_histogram(binwidth = bin_width, color = "black", alpha = 0.7, position = "identity") +
  labs(title = "P450 Gene Counts Across Lutzomyia longipalpis Genome", x = "Genome Position", y = "Count") +
  facet_wrap(~ seq_id, scales = "free_x", ncol = length(unique(genes$seq_id))) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

p450_localization


#facet_wrap(~ seq_id, scales = "free") +
#ggsave(plot=p450_localization, filename = "Ppapatasi_P450s_genome_localization.svg")
ggsave(plot=p450_localization, filename = "Llongipalpis_P450s_genome_localization2.svg")
