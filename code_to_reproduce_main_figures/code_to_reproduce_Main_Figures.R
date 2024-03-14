
source("~/bin/OliveTrees/OliveTree/R/OliveTree.R")
source("~/bin/OliveTrees/OliveTree/R/plotter.R")


#====================== FIGURE 1: Sand fly species tree and stacked barplots with orthogroup representation ===================#

# Figure 1A: Sand fly species tree
species_tree <- read_tree("Supplementary_File_4.Sandfly_species_tree.nwk")

species_tree@phylo$tip.label <- lapply(species_tree@phylo$tip.label, function(x) {
  gsub(".final_genesets|.geneset", "", x)
}) %>% unlist()

Figure1A <- visualize_tree(species_tree,
                           tiplabels = TRUE,
                           tip_label_size = 4,
                           bootstrap_numbers = FALSE,
                           bootstrap_circles = TRUE,
                           bootstrap_number_nudge_x = 0,
                           bootstrap_number_nudge_y = 0.2,
                           save = FALSE
                           )

ggsave(plot = Figure1A,
       filename = "Figure_1A.svg",
       dpi = 600)

# Figure 1B: Stacked barplots with orthogroup representation - Input is the orthology_results_for_R.txt produced by classify_orthogroups.pl

orthology_counts <- read.delim("orthology_results_for_R.txt", 
                               header = T
                               )

order_vector <- c(118:126, 136:144, 163:171, 73:81, 
                  82:90, 91:99, 127:135, 154:162,
                  109:117, 100:108, 145:153)

orthology_counts <- orthology_counts[order_vector,]

species_order <- unique(orthology_counts$Species)
orthology_counts$Species <- factor(orthology_counts$Species, 
                                   levels = unique_order)  # Set custom levels for the factor

stack_order <- rev(unique(orthology_counts$Type))[c(1,4:7,3,2,8,9)]
orthology_counts$Type <- factor(orthology_counts$Type, levels = stack_order)

Figure1B <- ggplot(orthology_counts, aes(orthology_counts =Species, y = Number_of_Genes, fill = Type)) + 
            geom_bar(position = "stack", stat = "identity") +
            scale_fill_manual(values = rev(c("darkblue", "skyblue2", "lightblue2", "orange1", "brown3", "salmon", "yellow2", "green4", "grey"))) +
            coord_flip() + theme_few() + scale_color_manual(values="black")

ggsave(plot = Figure1B,
       filename = "Figure_1B.svg",
       dpi = 600)


#====================== FIGURE 2: Total P450 phylogeny and P450 gene counts across species ===================#

# Figure 2A: P450 phylogeny with highlighted clans

P450_tree <- read_tree("Supplementary_File_6.Sandfly_1275_P450s_plus_Agambiae.nwk", 
                       bootstrap_support = TRUE)

P450_tree@phylo$tip.label <- lapply(P450_tree@phylo$tip.label, function(x) {
  gsub("^g", "longipalpis_g", x)
}) %>% unlist()

P450_tree@phylo$tip.label <- lapply(P450_tree@phylo$tip.label, function(x) {
  gsub("^Phlebotomus_papatasi", "papatasi", x)
}) %>% unlist()


# Use node_ids to identify the labels of parent nodes of the four P450 clans
node_ids(P450_tree, 
         form = "rectangular", 
         pattern = "Agam"
         )

# Nodes to highlight with corresponding CYP clan labels as names
highlight_nodes = c("MITO" = 1702,
                    "CYP2" = 1830,
                    "CYP3" = 1924,
                    "CYP4" = 1378)

# Generate the highlighted tree with colors corresponding to the different clans
Figure2A <- highlight_tree(P450_tree,
                           highlight_nodes = highlight_nodes, 
                           colors = c("CYP4"= "orange2",
                                      "CYP3" = "green4",
                                      "MITO" = "lightblue2",
                                      "CYP2" = "gold"
                                     )
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


Figure_S4 <- visualize_tree(tree = P450_tree,
                            form = "circular", 
                            tiplabels = TRUE,
                            pattern = names(group_colors),
                            tip_label_colors = group_colors, 
                            bootstrap_numbers = FALSE,
                            bootstrap_circles = TRUE,
                            )

# Figure 2B: Heatmap of P450 gene counts per CYP subfamily across the 11 sand fly species

source("complex_heatmap.R")

Figure2B <- complex_heatmap("P450_gene_counts.tsv", 
                            reordered_rows = c(2:7,1,8:58),
                            reordered_cols = c(12,5,3,7,2,1,4,6,8,9,10,11),
                            color_palette = c("azure2","dodgerblue3","dodgerblue4"), 
                            head_annotation =  c(MITO = "lightblue", CYP2 = "gold", CYP3 = "green4", CYP4 = "orange2"),
                            row_annotation = TRUE, 
                            gaps_row = c(8,10), 
                            gaps_col = c(7,15,44)
                            )
ggsave(plot = Figure2B,
       filename = "Figure2B.svg",
       dpi = 600)


#====================== FIGURE 3: Independent CYP expansions and P450 genomic clusters ===================#

# Figure 3A: Independent CYP expansions

CYP6ACJ <- extract_subtree(P450_tree, 
                           tip1 = "arabicus_TRINITY_DN1067_c1_g4_i1_p1_501", 
                           tip2 = "schwetzi_TRINITY_DN7802_c0_g1_i5_p1_511")

CYP9JR <- extract_subtree(P450_tree, 
                          tip1 = "schwetzi_NODE_11474_length_1711_cov_34_263736_g5680_i0_p1_529", 
                          tip2 = "longipalpis_g10133")

# Color and shape mappings for each sand fly species
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
                  schwetzi = "triangle star",
                  sergenti = "circle",
                  tobbi = "circle")


CYP6ACJ_tree <- visualize_tree(tree = CYP6ACJ, 
                               color = group_colors,
                               shape = group_shapes,
                               bootstrap_circles = TRUE,
                               bootstrap_numbers = FALSE,
                               tiplabels = FALSE, 
                               tip_label_size = 4,
                               branch_length = FALSE)

CYP9JR_tree <- visualize_tree(tree = CYP9JR, 
                              color = group_colors,
                              shape = group_shapes,
                              bootstrap_circles = TRUE, 
                              bootstrap_numbers = FALSE,
                              tiplabels = FALSE, 
                              tip_label_size = 4,
                              branch_length = FALSE)

Figure3A <- multipanel(CYP6ACJ_tree, 
                       CYP9JR_tree, 
                       horizontal = TRUE)


# Figure 3B: P450 clusters in Ph. papatasi and Lu. longipalpis genomes

source("gggenomes_lib.R")

CYP6ACJ_cluster <- visualize_clusters(gff3_file = "CYP6ACJ.clusters.gff3")

CYP9JR_cluster <- visualize_clusters(gff3_file = "CYP9JR.clusters.gff3")

Figure3B <- multipanel(CYP6ACJ_cluster,
                       CYP9JR_cluster,
                       horizontal = FALSE)

# Combine Figure3A and Figure3B into a single figure 
Figure3 <- multipanel(Figure3A,
                      Figure3B,
                      horizontal = FALSE)

ggsave(plot = Figure3,
       filename = "Figure_3.svg",
       dpi = 600)

#=================================== FIGURE 4: CYP4G347 is probably the functional ortholog of CYP4G17 =====================================#

# Figure 4A: Sand flies have two CYP4G17 orthologs
CYP4G17 <- extract_subtree(P450_tree, 
                           tip1 = "Agam_4G17_clan4_562", 
                           tip2 = "papatasi_CYP4G347_3_13858k")

node_labels = c("CYP4G346" = 38, 
                "CYP4G347" = 28)

Figure4A <- visualize_tree(tree = CYP4G17, 
                           color = group_colors,
                           shape = group_shapes,
                           bootstrap_circles = TRUE, 
                           bootstrap_numbers = FALSE,
                           tiplabels = FALSE, 
                           tip_label_size = 4,
                           branch_length = FALSE,
                           reference = "Agam")
ggsave(plot = Figure4A,
       filename = "Figure4A.svg",
       dpi = 600 )


# Figure 4B: Violin plots with expression of P450 gene expression across the 11 sand fly species

source("violin_plots.R")

CYP4G346 <- extract_subtree(CYP4G17, 
                            tip1 = "duboscqi_TRINITY_DN414",
                            tip2 = "argentipes_TRINITY_DN13789"
                            )

CYP4G347 <- extract_subtree(CYP4G17, 
                            tip1 = "tobbi_NODE_14046",
                            tip2 = "sergenti_NODE_11431"
                            )

CYP4G346_orthologs <- lapply (CYP4G346@phylo$tip.label, function (x) { x = "red"} )

CYP4G346_orthologs <- setNames(rep("red", length(CYP4G346@phylo$tip.label)), CYP4G346@phylo$tip.label)
CYP4G347_orthologs <- setNames(rep("green", length(CYP4G347@phylo$tip.label)), CYP4G347@phylo$tip.label)

sgenes <- c(CYP4G346_orthologs, CYP4G347_orthologs)

Figure4B <- expression_violin_plot("all.P450s.tpm3", 
                                   title = "CYP4G346 is the most highly expressed P450 gene", 
                                   selected_genes = sgenes, 
                                   mode = "mean",
                                   parent_color = "grey", 
                                   )

Figure4 <- multipanel(Figure4A, 
                      Figure4B,
                      horizontal = FALSE
                      )

#======================================================== FIGURES 5: GST variation is mostly located on Lutzomyia-specific GSTD and GSTX expansions ========================================================#

# Figure 5A: Heatmap of GST gene counts per species
Figure5A <- complex_heatmap("GST_gene_counts.tsv", 
                            color_palette = c("azure2","dodgerblue3","dodgerblue4"), 
                            row_annotation = TRUE, 
                            gaps_row = c(8,10)
                            )

ggsave(plot = Figure5A, 
       filename = "Figure_5A_GST_counts_heatmap.svg", 
       dpi = 600)


# Figure 5B: Independent GSTD and GSTX expansions
GST_tree <- read_tree("Sandfly_309_GSTs_plus_Agambiae.nwk")
GSTD <- extract_subtree(GST_tree, "AGAP004165_GSTd2", "schwetzi_TRINITY_DN300_c0_g1_i18_p1")
GSTX <- extract_subtree(GST_tree, "schwetzi_TRINITY_DN5253_c0_g2_i4_p1", "schwetzi_NODE_7136_length_2426_cov_49_318317_g936_i5_p1")

GSTD_tree <- visualize_tree(tree = GSTD, 
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_legend = TRUE,
                            bootstrap_numbers = FALSE,
                            tiplabels = FALSE, 
                            tip_label_size = 6,
                            branch_length = FALSE,
                            reference1 = "AGAP",
                            save = FALSE
                            )

GSTX_tree <- visualize_tree(tree = GSTX, 
                            color = group_colors,
                            shape = group_shapes,
                            bootstrap_circles = TRUE, 
                            bootstrap_legend = TRUE,
                            bootstrap_numbers = FALSE,
                            tiplabels = FALSE, 
                            tip_label_size = 6,
                            branch_length = FALSE,
                            reference1 = "AGAP",
                            save = FALSE
                            )

Figure5B <- multipanel(GSTD_tree, 
                       GSTX_tree, 
                       horizontal = TRUE)

ggsave(plot = Figure5B, 
       filename = "Figure_5B_GST_independent_expansions.svg", 
       dpi = 600)


Figure5 <- multipanel(Figure5A,
                      Figure5B,
                      horizontal = FALSE
                     )


Figure_S7 <- visualize_tree(tree = GST_tree,
                            form = "circular", 
                            tiplabels = TRUE,
                            pattern = names(group_colors),
                            tip_label_colors = group_colors, 
                            bootstrap_numbers = FALSE,
                            bootstrap_circles = TRUE,
)


#======================================================== FIGURES 6A, 6B: Sand fly orthologs for Ace1 and ABCB Mdr Full Transporters ==========================================================#

CCE_tree <- read_tree("/home/jason/Documents/sandflies/submission/Main_Figures/all_genes.fasta.trimmed.aln.phy.raxml.support.tree", 
                      bootstrap_support = TRUE)

# Extract subtree anchored by the An. gambiae Ace1(AGAP001356) and Ace2(AGAP00466)
ACE <- extract_subtree(CCE_tree, 
                       tip1 = "AGAP001356", 
                       tip2 = "AGAP000466")

# Identify node numbers corresponding to Ace1 and Ace2
node_ids(ACE, pattern = "AGAP")

clade_labs <- c("Ace1" = 37, 
                "Ace2" = 26)

# Visualize the subtree with the Ace1 and Ace2 orthologs, with color/shape species mappings and bootstrap colored circles and AGAP-matching IDs are reference
ACE_tree <- visualize_tree(tree = ACE, 
                           color = group_colors, 
                           shape = group_shapes, 
                           bootstrap_circles = TRUE, 
                           bootstrap_numbers = FALSE,
                           branch_length = FALSE,
                           tip_label_size = 3,
                           clades = clade_labs,
                           save = FALSE,
                           reference = "AGAP")

# Visualize ABCB FT tree with color/shape species mappings and bootstrap colored circles
ABCB <- read_tree("abcbs.nwk")

node_ids(ABCB, pattern = "Dmel")

clade_labs <- c("Mdr49" = 54, 
                "Mdr65" = 67,
                "Mdr50" = 42)

# Visualize the subtree with the Ace1 and Ace2 orthologs, with color/shape species mappings and bootstrap colored circles and AGAP-matching IDs are reference
ABCB_tree <- visualize_tree(tree = ABCB, 
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
                            reference1 = "Dmel",
                            reference2 = "Agam"
                            )

# Combine Figures horizontally
Figure6AB <- multipanel(ACE_tree,
                        ABCB_tree, 
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
L995 <- visualize_msa(aln_filename = "VGSC_orthologs.aln", range = c(1021, 1030))
N1570 <- visualize_msa(aln_filename = "VGSC_orthologs.aln", range = c(1604, 1614))

Figure6C <- multipanel(L995, 
                       N1570, 
                       horizontal = TRUE)

Figure6 <- multipanel(Figure6AB, 
                      Figure6C, 
                      horizontal = FALSE)

ggsave(plot = Figure6, 
       filename = "Figure_6.svg", 
       dpi = 600)



