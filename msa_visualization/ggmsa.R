library("ggmsa")
library("ggplot2")

# Define a function to visualize MSA using ggmsa
visualize_msa <- function(aln_filename, range, output_name="MSA") {
  
  ggmsa(aln_filename, range[1], range[2], seq_name = TRUE, char_width = 0.5) +
  geom_seqlogo(color = "Chemistry_AA") +
  geom_msaBar()
  
}

# Call the function with the specified alignment file ("ace1.aln") and range (c(250, 300))
l995 <- visualize_msa(aln_filename = "VGSC_orthologs.aln", range = c(1021, 1030))

n1570 <- visualize_msa(aln_filename = "VGSC_orthologs.all.aln", range = c(1604, 1614))

ace1 <- visualize_msa(aln_filename = "ace1.all.fasta.aln", range=c(280,310))

ggsave(plot=l995, file="L995.svg")
ggsave(plot=n1570, file="N1570.svg")
ggsave(plot=ace1, file="ACE1.svg")
