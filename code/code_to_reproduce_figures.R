
source("~/bin/OliveTrees/OliveTree/R/OliveTree.R")
source("~/bin/OliveTrees/OliveTree/R/plotter.R")

setwd("~/Documents/sandflies/submission/supplemental_information/Supplementary_Files")

#====================== FIGURE 1: Sand fly species tree and stacked barplots with orthogroup representation ===================#

# Figure_1A: Sand fly species tree
species_tree <- read_tree("Supplementary_File_4.Sandfly_species_tree.nwk")

species_tree@phylo$tip.label <- lapply(species_tree@phylo$tip.label, function(x) {
  gsub(".final_genesets|.geneset", "", x)
}) %>% unlist()

visualize_tree( tree = species_tree,
                tiplabels = TRUE,
                tip_label_size = 4,
                bootstrap_numbers = TRUE,
                bootstrap_circles = FALSE,
                bootstrap_number_nudge_x = 0,
                bootstrap_number_nudge_y = 0.4,
                save = TRUE,
                output = "Figure_1A.svg"
               )

# Figure_1B: Stacked barplots with orthogroup representation - Input is the orthology_results_for_R.txt produced by classify_orthogroups.pl
library(ggthemes)

orthology_counts <- read.delim("Supplementary_File_5.Orthology_results_for_ggplot.txt", header = T)

order_vector <- c(118:126, 136:144, 163:171, 73:81, 
                  82:90, 91:99, 127:135, 154:162,
                  109:117, 100:108, 145:153)

orthology_counts <- orthology_counts[order_vector,]

species_order <- unique(orthology_counts$Species)
orthology_counts$Species <- factor(orthology_counts$Species, 
                                   levels = species_order)  # Set custom levels for the factor

stack_order <- rev(unique(orthology_counts$Type))[c(1,4:7,3,2,8,9)]
orthology_counts$Type <- factor(orthology_counts$Type, levels = stack_order)

Figure_1B <- ggplot(orthology_counts, aes(x = Species, y = Number_of_Genes, fill = Type)) + 
             geom_bar(position = "stack", stat = "identity") +
             scale_fill_manual(values = rev(c("darkblue", "skyblue2", "lightblue2", "orange1", "brown3", "salmon", "yellow2", "green4", "grey"))) +
             coord_flip() + theme_few() + scale_color_manual(values="black")

ggsave(plot = Figure_1B,
       filename = "Figure_1B.svg",
       dpi = 600)


#====================== FIGURE 2: Total P450 phylogeny and P450 gene counts across species ===================#

# Figure_2A: P450 phylogeny with highlighted CYP clans
P450_tree <- read_tree("Supplementary_File_7.Sandfly_1275_P450s_plus_Agambiae.nwk", bootstrap_support = TRUE)

P450_tree@phylo$tip.label <- unlist(lapply(P450_tree@phylo$tip.label, function(x) {
  gsub("^g", "longipalpis_g", x)
}))

P450_tree@phylo$tip.label <- unlist(lapply(P450_tree@phylo$tip.label, function(x) {
  gsub("^Phlebotomus_papatasi", "papatasi", x)
}))


# Use print_internal_nodes to identify the labels of parent nodes of the four P450 clans
print_internal_nodes(P450_tree, 
                     form = "circular", 
                     references = "Agam")

# Nodes to highlight with corresponding CYP clan labels as names
highlight_nodes = c("MITO" = 1702,
                    "CYP2" = 1830,
                    "CYP3" = 1924,
                    "CYP4" = 1378)

# Generate the highlighted tree with colors corresponding to the different clans
Figure_2A <- highlight_tree(P450_tree,
                            form = "rectangular",
                            highlight_nodes = highlight_nodes, 
                            colors = c("CYP4"= "orange2",
                                       "CYP3" = "green4",
                                       "MITO" = "lightblue2",
                                       "CYP2" = "gold")
                            )

# Color mappings for each sand fly species
group_colors <- c(arabicus = "purple",
                  argentipes =  "cyan",
                  duboscqi = "salmon2",
                  longipalpis = "darkgreen",
                  migonei =  "green",
                  orientalis =  "gold",
                  perniciosus =  "lightgreen",
                  papatasi =  "darkred",
                  schwetzi  = "blue",
                  sergenti =  "orange",
                  tobbi =  "red")

group_shapes <- c(arabicus = "circle",
                  argentipes = "circle",
                  duboscqi = "circle",
                  longipalpis = "square",
                  migonei = "square",
                  orientalis = "circle",
                  papatasi = "circle",
                  perniciosus = "circle",
                  schwetzi = "hexagonal star",
                  sergenti = "circle",
                  tobbi = "circle")


# Figure_S4: Visualize P450 tree with reference tip labels (CYP names) from Ph. papatasi and An. gambiae
visualize_tree(tree = P450_tree,
               form = "rectangular", 
               #color = group_colors,
               #shape = group_shapes,
               #references = c("papatasi", "Agam"),
               tip_label_size = 2,
               tip_shape_size = 0.8,
               #tip_label_colors = group_colors, 
               tiplabels = TRUE,
               pattern_id = "Agam",
               bootstrap_numbers = FALSE,
               bootstrap_circles = TRUE,
               bootstrap_circle_size = 0.7
               save = TRUE,
               output = "Figure_S4_P450_phylogeny_9_transcriptomes_2_genomes.svg"
               )


# Figure_2B: Heatmap of P450 gene counts per CYP subfamily across the 11 sand fly species
source("~/bin/SandFlyComparativeGenomics/code/libraries/complex_heatmap.R")

Figure_2B <- complex_heatmap("P450_gene_counts.tsv", 
                             reordered_rows = c(2:7,1,8:58),
                             reordered_cols = c(12,5,3,7,2,1,4,6,8,9,10,11),
                             color_palette = c("azure2","dodgerblue3","dodgerblue4"), 
                             head_annotation =  c(MITO = "lightblue", CYP2 = "gold", CYP3 = "green4", CYP4 = "orange2"),
                             row_annotation = TRUE, 
                             gaps_row = c(8,10), 
                             gaps_col = c(7,15,44)
                             )

# Combine Figure_2A and Figure_2B into a single figure vertically
Figure_2 <- multipanel(Figure_2A,
                       Figure_2B,
                       horizontal = FALSE)

ggsave(plot = Figure_2B,
       filename = "Figure_2B.svg",
       dpi = 600)

#====================== FIGURE 3: Independent CYP expansions and P450 genomic clusters ===================#

# Figure_3A: Independent CYP expansions
CYP6ACJ <- extract_subtree(P450_tree, 
                           tip1 = "arabicus_TRINITY_DN1067_c1_g4_i1_p1_501", 
                           tip2 = "schwetzi_TRINITY_DN7802_c0_g1_i5_p1_511")

CYP9JR <- extract_subtree(P450_tree, 
                          tip1 = "schwetzi_NODE_11474_length_1711_cov_34_263736_g5680_i0_p1_529", 
                          tip2 = "longipalpis_g10133")

CYP6ACJ_tree <- visualize_tree(tree = CYP6ACJ, 
                               color = group_colors,
                               shape = group_shapes,
                               bootstrap_numbers = FALSE,
                               tiplabels = FALSE, 
                               tip_label_size = 4,
                               branch_length = TRUE)

CYP9JR_tree <- visualize_tree(tree = CYP9JR, 
                              color = group_colors,
                              shape = group_shapes,
                              bootstrap_circles = TRUE, 
                              bootstrap_numbers = FALSE,
                              tiplabels = FALSE, 
                              tip_label_size = 4,
                              branch_length = FALSE)

Figure_3A <- multipanel(CYP6ACJ_tree, 
                        CYP9JR_tree, 
                        horizontal = TRUE)

# Figure_3B: P450 clusters in Ph. papatasi and Lu. longipalpis genomes
source("~/bin/SandFlyComparativeGenomics/code/libraries/gggenomes_lib.R")

CYP6ACJ_cluster <- visualize_clusters(gff3_file = "CYP6ACJ.clusters.gff3")
CYP9JR_cluster <- visualize_clusters(gff3_file = "CYP9JR.clusters.gff3")

Figure_3B <- multipanel(CYP6ACJ_cluster,
                        CYP9JR_cluster,
                        horizontal = FALSE)

# Combine Figure_3A and Figure_3B into a single figure 
Figure_3 <- multipanel(Figure_3A,
                       Figure_3B,
                       horizontal = FALSE)

ggsave(plot = Figure_3,
       filename = "Figure_3.svg",
       dpi = 600)


# Figure_S5: Compute gene counts across Ph. papatasi and Lu. longipalpis genomes
source("~/bin/SandFlyComparativeGenomics/code/libraries/count_genes_sliding_window.R")

Figure_S5 <- count_genes_sliding_window(coordinate_file = "P450_genes_in_sandfly_genomes.bed",
                                        sliding_window_size = 50000)
ggsave(plot = Figure_S5,
       filename = "Figure_S5.svg",
       dpi = 600)

#=================================== FIGURE 4: CYP4G347 is probably the functional ortholog of CYP4G17 =====================================#

# Figure_4A: Sand flies have two CYP4G17 orthologs
CYP4G17_tree <- extract_subtree(P450_tree, 
                                tip1 = "Agam_4G17_clan4_562", 
                                tip2 = "papatasi_CYP4G347_3_13858k")

node_labels = c("CYP4G346" = 38, 
                "CYP4G347" = 28)

Figure_4A <- visualize_tree(tree = CYP4G17_tree, 
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_numbers = FALSE,
                            tiplabels = FALSE, 
                            tip_label_size = 4,
                            branch_length = FALSE,
                            references = "Agam")
ggsave(plot = Figure_4A,
       filename = "Figure_4A.svg",
       dpi = 600 )


# Figure_4B: Violin plots with expression of P450 gene expression across the 11 sand fly species
source("~/bin/SandFlyComparativeGenomics/code/libraries/violin_plots.R")

CYP4G346 <- extract_subtree(CYP4G17, 
                            tip1 = "duboscqi_TRINITY_DN414",
                            tip2 = "argentipes_TRINITY_DN13789")

CYP4G347 <- extract_subtree(CYP4G17, 
                            tip1 = "tobbi_NODE_14046",
                            tip2 = "sergenti_NODE_11431")

CYP4G346_orthologs <- setNames(rep("red", length(CYP4G346@phylo$tip.label)), CYP4G346@phylo$tip.label)
CYP4G347_orthologs <- setNames(rep("green", length(CYP4G347@phylo$tip.label)), CYP4G347@phylo$tip.label)

sgenes <- c(CYP4G346_orthologs, CYP4G347_orthologs)

Figure_4B <- expression_violin_plot("all.P450s.tpm", 
                                    title = "CYP4G346 is the most highly expressed P450 gene", 
                                    selected_genes = sgenes, 
                                    mode = "mean",
                                    parent_color = "grey")

Figure_4 <- multipanel(Figure_4A, 
                       Figure_4B,
                       horizontal = FALSE)


#======================================================== FIGURES 5: GST variation is mostly located on Lutzomyia-specific GSTD and GSTX expansions ========================================================#

# Figure_5A: Heatmap of GST gene counts per species
Figure_5A <- complex_heatmap(input_file = "GST_gene_counts.tsv", 
                             color_palette = c("azure2","dodgerblue3","dodgerblue4"),
                             reordered_cols = c(5,3,7,2,1,4,6,8,9,10,11,12),
                             right_annotation = TRUE, 
                             right_annotation_title = "Total genes",
                             left_annotation = TRUE,
                             left_annotation_title = "Species",
                             gaps_row = c(8,10),
                             fontsize = 10,
                             legend = TRUE,
                             legend_title = "Gene counts"
                            )

ggsave(plot = Figure_5A, 
       filename = "Figure_5A_GST_counts_heatmap.svg", 
       dpi = 600)


# Figure_5B: Independent GSTD and GSTX expansions
GST_tree <- read_tree("Supplementary_File_9.Sandfly_276_GSTs_plus_Agambiae_and_Hassan_characterized_GSTs.nwk")

GSTD <- extract_subtree(GST_tree, "AGAP004165_GSTd2", "schwetzi_TRINITY_DN300_c0_g1_i18_p1")
GSTX <- extract_subtree(GST_tree, "migonei_TRINITY_DN6514_c0_g1_i16_p1", "perniciosus_TRINITY_DN1375_c0_g2_i1_p1")

GSTD_tree <- visualize_tree(tree = GSTD, 
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_legend = TRUE,
                            bootstrap_numbers = FALSE,
                            tiplabels = FALSE, 
                            tip_label_size = 3,
                            branch_length = FALSE,
                            references = "AGAP",
                            save = FALSE)

GSTX_tree <- visualize_tree(tree = GSTX, 
                            flip_nodes = TRUE,
                            node1 = 46,
                            node2 = 68,
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_legend = TRUE,
                            bootstrap_numbers = FALSE,
                            tiplabels = FALSE, 
                            tip_label_size = 3,
                            branch_length = FALSE,
                            references = "Hassan",
                            save = FALSE)

Figure_5B <- multipanel(GSTD_tree, 
                        GSTX_tree, 
                        horizontal = TRUE)

ggsave(plot = Figure_5B, 
       filename = "Figure_5B_GST_independent_expansions.svg", 
       dpi = 600)

Figure_5 <- multipanel(Figure_5A,
                       Figure_5B,
                       horizontal = FALSE)

# Visualize GST tree
group_colors <- c(group_colors, 
                  "AGAP" = "black",
                  "Hassan" = "black")

Figure_S10 <- visualize_tree(tree = GST_tree,
                            form = "circular", 
                            color = group_colors,
                            references = c("AGAP", "Hassan"),
                            tip_shape_size = 0.8,
                            tip_label_size = 1.5,
                            bootstrap_numbers = FALSE,
                            bootstrap_circles = TRUE,
                            bootstrap_circle_size = 0.8,
                            bootstrap_legend = TRUE,
                            save = TRUE,
                            output = "Figure_S10_276_GSTs_phylogeny_plus_Agambiae_and_Hassan_characterized_GSTs.svg")


# Read and Visualize UGT tree
UGT_tree <- read_tree("Supplementary_File_11.Sandfly_214_UGTs_plus_Agambiae.nwk")

Figure_S9 <- visualize_tree(tree = UGT_tree,
                            form = "circular", 
                            tiplabels = TRUE,
                            tip_label_colors = group_colors,
                            pattern_id = names(group_colors),
                            tip_shape_size = 2.5,
                            mappings_legend = FALSE,
                            tip_label_size = 2,
                            bootstrap_numbers = FALSE,
                            bootstrap_circles = TRUE,
                            bootstrap_circle_size = 0.8,
                            bootstrap_legend = FALSE,
                            save = TRUE,
                            output = "Figure_S9_UGTs_phylogeny.svg")


# Read and Visualize CCE tree
CCE_tree <- read_tree("Supplementary_File_13.Sandfly_379_CCEs_plus_Agambiae.nwk", bootstrap_support = TRUE)

Figure_S11 <- visualize_tree(CCE_tree, 
                             form = "circular",
                             tiplabels = TRUE,
                             tip_shape_size = 0.8,
                             references = "AGAP",
                             tip_label_size = 1,
                             color = group_colors,
                             shape = group_shapes,
                             bootstrap_circles = TRUE,
                             bootstrap_numbers = FALSE,
                             bootstrap_circle_size = 0.8,
                             save = TRUE,
                             output = "Figure_S11_CCE_phylogeny.svg")


# Extract subtree anchored by the An. gambiae Ace1 and Ace2 genes
ACE_tree <- extract_subtree(CCE_tree, 
                            tip1 = "AGAP001356_ACE1", 
                            tip2 = "AGAP000466_ACE2")

# Identify node numbers corresponding to Ace1 and Ace2
print_internal_nodes(ACE_tree, pattern = "AGAP")

clade_labs <- c("Ace1" = 37, 
                "Ace2" = 26)


# Visualize the subtree with the Ace1 and Ace2 orthologs, with color/shape species mappings and bootstrap colored circles and AGAP-matching IDs are reference
Figure_6A <- visualize_tree(tree = ACE_tree, 
                            color = group_colors, 
                            shape = group_shapes, 
                            bootstrap_circles = TRUE, 
                            bootstrap_numbers = FALSE,
                            branch_length = FALSE,
                            tip_label_size = 3,
                            clades = clade_labs,
                            save = FALSE,
                            references = "AGAP")


# Visualize ABCB FT tree with color/shape species mappings and bootstrap colored circles
ABC_tree <- read_tree("Supplementary_File_15.Sandfly_561_ABC_transporters_plus_Agambiae.nwk")

Figure_S13 <- visualize_tree(tree = ABC_tree,
                             form = "circular", 
                             tiplabels = TRUE,
                             color = group_colors,
                             shape = group_shapes,
                             tip_shape_size = 0.8,
                             references = "AGAP",
                             mappings_legend = FALSE,
                             tip_label_size = 1.5,
                             bootstrap_numbers = FALSE,
                             bootstrap_circles = TRUE,
                             bootstrap_circle_size = 0.8,
                             bootstrap_legend = FALSE,
                             save = TRUE,
                             output = "Figure_S13_ABC_phylogeny.svg")

# Visualize ABCB FT tree with color/shape species mappings and bootstrap colored circles
ABCB_tree <- read_tree("Supplementary_File_16.Sandfly_ABCB_transporters_plus_Agambiae_Dmel.nwk", bootstrap_support = TRUE)

print_internal_nodes(ABCB_tree, pattern = "Dmel")

clade_labs <- c("Mdr49" = 54, 
                "Mdr65" = 67,
                "Mdr50" = 42)

Figure_6B <- visualize_tree(tree = ABCB_tree, 
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_legend = TRUE,
                            bootstrap_numbers = FALSE,
                            bootstrap_number_nudge_x = 0.05,
                            bootstrap_number_nudge_y = 0.3,
                            tiplabels = FALSE, 
                            tip_label_size = 4,
                            branch_length = FALSE,
                            clades = clade_labs,
                            labeldist = 0.5,
                            references = c("Dmel", "Agam")
                            )

# Combine Figures 6A and 6B horizontally
Figure_6AB <- multipanel(ACE_tree,
                         Figure_6B, 
                         horizontal = "T")

#======================================================== FIGURE 6C: Multiple Sequence Alignment of VGSC orthologs in critical, resistance-conferring regions ==========================================================#
library("ggmsa")

# Define a function to visualize MSA using ggmsa
visualize_msa <- function(aln_filename, range, output_name = "MSA") {
  
  ggmsa(aln_filename, range[1], range[2], seq_name = TRUE, char_width = 0.5) +
    geom_seqlogo(color = "Chemistry_AA") +
    geom_msaBar()
  
}

# Call the function with the specified alignment file("ace1.aln") and range(c(250, 300))
L995 <- visualize_msa(aln_filename = "Supplementary_File_18.Sandfly_11_VGSC_orthologs_plus_Agambiae.aln", range = c(1021, 1030))
N1570 <- visualize_msa(aln_filename = "Supplementary_File_18.Sandfly_11_VGSC_orthologs_plus_Agambiae.aln", range = c(1604, 1614))

Figure_6C <- multipanel(L995, 
                        N1570, 
                        horizontal = TRUE)

Figure_6 <- multipanel(Figure_6AB, 
                       Figure_6C, 
                       horizontal = FALSE)

ggsave(plot = Figure_6, 
       filename = "Figure_6.svg", 
       dpi = 600)
