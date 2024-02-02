#Library of functions for advanced tree manipulation and visualization using ggtree, ape, phytools, and other related tools.

# Function to check if a package is installed, and if not, install it.
package_install <- function(package_name) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  } else {
    print(sprintf("%s %s", package_name, "is not installed. Installing it!"))
    is_available <- BiocManager::available(package_name)
    
    if (any(is_available == "TRUE")) {
      BiocManager::install(package_name, dependencies = TRUE, update = TRUE)
      } else {
        install.packages(package_name, dependencies = TRUE, ask = FALSE, reinstall = TRUE)
    }
  }
}

find_previous_versions <- function () {
  tmp = as.data.frame(installed.packages()) 
  max_version = max(as.numeric(substr(tmp$Built, 1,1)))
  tmp = tmp[as.numeric(substr(tmp$Built, 1,1)) < max_version,]
  lapply(tmp$Package, remove.packages)
  lapply(tmp$Package, function(x) install.packages(x, dependencies = TRUE))
}
BiocManager::install("devtools")

library(BiocManager)
BiocManager::install("ggmsa")


find_previous_versions()

# Load required packages or install them if necessary
dependencies <- c("ape", "phytools", "treeio", "TreeTools", "ggstar", 
                  "ggtree", "ggplot2", "dplyr", "stringi", "stringr")

for (pkg in dependencies) {
  package_install(pkg)
}

install.packages(names(dependencies_list), dependencies = TRUE)


#==================================== TREE MANIPULATION FUNCTIONS ====================================#

# Function to read and preprocess a tree
read_tree <- function(nwk) {
  tree <- ape::read.tree(nwk)
  return(Preorder(tree))
}

# Function to print a tree with node IDs
node_ids <- function(tree, reference = NULL) {
  # Create the base tree plot with node labels
  tree_plot <- ggtree(tree, layout = "rectangular") + 
    geom_nodelab(aes(label = node), hjust = -0.1, color = "red")
  
  if (is.null(reference)) {
    return(tree_plot) 
  } else {
    # Filter the tips based on the reference string
    tips_to_label <- tree$tip.label[grep(reference, tree$tip.label)]
    
    # Create the plot with tip labels for the specified tips
    labeled_tree <- tree_plot + geom_tiplab2(
      aes(subset = label %in% tips_to_label, label = label, color = "red"), 
      geom = "label"
    )
    
    return(labeled_tree)
  }
}

# Function to collapse nodes based on bootstrap support, returns a "phylo" object
bootstrap_collapse <- function(tree, cutoff) { 
  return(as.polytomy(tree, feature = 'node.label', fun = function(x) as.numeric(x) < cutoff))
}

# Function to flip based on descendant nodes (internal nodes or leaves if node is terminal)
flip_node <- function(tree, node1, node2) {
  return(as.phylo(flip(ggtree(tree), node1, node2)))
}

# Function to group all descendant branches of node(s)
group_descendants <- function(tree, node1, node2 = "", node3 = "", node4 = "") {
  return(groupClade(tree, .node = c(node1, node2, node3, node4)))
}

# Function to extract a subtree by finding the MRCA of two anchor nodes while preserving branch lengths and bootstrap values 
extract_subtree <- function(t, tip1, tip2, bl = TRUE) {
  if (!(bl == "TRUE" || bl == "T")) {
    t$edge.length <- NULL
  } 
  return(ape::extract.clade(t, node = phytools::findMRCA(t, c(tip1, tip2))))
}

#==================================== TREE VISUALIZATION FUNCTIONS ====================================#

# Function to include bootstrap values as colored circles
bootstrap_circles <- function(x) {
  bootstrap_colors <- c('(75,100]' = "black", '(50,75]' = "grey", '(0,50]' = "snow2")
  categories <- c('(75,100]', '(50,75]', '(0,50]')
  
  bs_tibble <- tibble(
    node = 1:Nnode(x@phylo) + Ntip(x@phylo),
    bootstrap = as.numeric(x@phylo$node.label),
    cut(bootstrap, c(0, 50, 75, 100))
  ) 
  
  return(data.frame(
    node = bs_tibble$node,
    bootstrap = as.numeric(bs_tibble$bootstrap),
    category = cut(bs_tibble$bootstrap, c(0, 50, 75, 100)),
    b_color = bootstrap_colors[match(cut(bs_tibble$bootstrap, c(0, 50, 75, 100)), categories)]
  )) 
}

# Function to highlight nodes on a tree
highlight_tree <- function(tree, highlight_nodes, colors = NULL, layout = "circular", name = NULL, ...) {
  # Check if highlight_nodes is a list or a vector
  if (is.list(highlight_nodes)) {
    # It's a list, create a data frame with random colors
    if (is.null(colors)) {
      # Generate random colors if colors are not provided
      highlight <- data.frame(
        Groups = names(highlight_nodes),  
        Label = highlight_nodes,
        Color = sample(colors(), length(highlight_nodes))
      )
    } else {
      # Use provided colors
      highlight <- data.frame(
        Groups = names(highlight_nodes),  
        Label = highlight_nodes,
        Color = colors
      )
    }
  } else if (is.vector(highlight_nodes)) {
    # It's a vector, create a data frame with provided colors or generate random colors
    if (is.null(colors)) {
      # Generate random colors if colors are not provided
      highlight <- data.frame(
        Groups = NULL,  # You can specify the groups if needed
        Label = highlight_nodes,
        Color = sample(colors(), length(highlight_nodes))
      )
    } else {
      # Use provided colors
      highlight <- data.frame(
        Groups = NULL,  # You can specify the groups if needed
        Label = highlight_nodes,
        Color = colors
      )
    }
  } else {
    stop("highlight_nodes must be either a list or a vector.")
  }
  
  # Create the ggtree plot
  if (any(grepl(".nwk|.tre", tree))) {
    x <- read_tree(tree)
    tree <- ggtree(x, layout = layout) + 
      geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) + 
      scale_starshape_identity() + scale_fill_identity()
  } else {
    tree <- ggtree(tree, layout = layout) + 
      geom_star(aes(subset = isTip, starshape = "circle"), fill = "lightgrey", size = 0.8) + 
      scale_starshape_identity() + scale_fill_identity()
  }
  
  # Highlight the nodes
  if (!is.null(highlight$Groups)) {
    treef <- tree + geom_hilight(
      data = highlight,
      aes(node = Label, fill = Groups),
      alpha = 0.3, extend = 0.10, linetype = 1, linewidth = 0.9, 
      colour = "black", show.legend = TRUE
    ) + scale_fill_manual(values = highlight$Color)
  } else {
    treef <- tree + geom_hilight(
      data = highlight,
      aes(node = Label, fill = highlight$Color),
      alpha = 0.3, extend = 0.10, linetype = 1, linewidth = 0.9, 
      colour = "black", show.legend = TRUE
    ) 
  }
  
  # Save the highlighted plot
  if (exists("treef")) {
    if (any(grepl(".nwk|.tre", tree))) {
      ggsave(plot = treef, sprintf("%s_highlighted.svg", sub(".nwk|.tre", "", tree)), dpi = 600)
      sprintf("Tree plotted and saved as %s_highlighted.svg", sub(".nwk|.tre", "", tree))
    } else {
      ggsave(plot = treef, "tree_plot_highlighted.svg", dpi = 600)
      print("Tree plotted and saved as tree_plot_highlighted.svg!")
    }
    return(treef)
  }
}

# Function to visualize a tree with a wide variety of options and customizable features
visualize_tree <- function(tree, color = NULL, shape = NULL, node1 = NULL, node2 = NULL, node3 = NULL, 
                           reference1 = NULL, reference2 = NULL, bootstrap_legend = TRUE, 
                           clades = NULL, labels = NULL, fontsize = 7, labeldist = 0.1, ...) {
  # Check if species matches reference and discard from the rest of tip labels - at least two references, adjust for more
  if (is.null(reference1) & is.null(reference2)) {
    ref_labs <- NULL
    ref_species <- NULL
    ref_idx <- NULL
    species <- sub("_.*", "", tree$tip.label)
  } else {  
    if (!is.null(reference2)) {
      reference <- c(reference1, reference2)
      ref_idx <- grep(paste(reference, collapse = "|"), tree$tip.label) 
    }
    if (is.null(reference2) & !is.null(reference1)) {
      reference <- reference1
      ref_idx <- grep(paste(reference), tree$tip.label) 
    }
    ref_labs <- tree$tip.label[ref_idx]
    ref_species <- data.frame(node = ref_idx, name = c(sub("_.*", "", ref_labs)), label = ref_labs)
    species <- sub("_.*", "", tree$tip.label[-ref_idx])
  }
  
  # Get index of identified species
  species_idx <- grep(paste(species, collapse = "|"), tree$tip.label)
  
  # Keep name for tip color and shape mapping per species
  unique_species <- sort(unique(species))
  
  # Create associative dataframes of tip color and shape 
  tip_colors_df <- data.frame(
    label = tree$tip.label[species_idx],
    species = species,
    colour = color[match(species, unique_species)]
  )
  
  tip_shapes_df <- data.frame(
    label = tree$tip.label[species_idx],
    species = species,
    shape = shape[match(species, unique_species)]
  )
  
  # Join tip colors and shapes in tree object based on tip label
  x <- data.frame()
  x <- full_join(tree, tip_colors_df, by = 'label')
  x <- full_join(x, tip_shapes_df, by = 'label')
  
  # Include bootstrap values as circles 
  bs_values <- bootstrap_circles(x)
  x <- full_join(x, bs_values, by = 'node')
  
  bootstrap_colors <- c('(75,100]' = "black", '(50,75]' = "grey", '(0,50]' = "snow2")
  
  if (bootstrap_legend == TRUE) {
    bootstrap_legend = 'legend'
  } else {
    bootstrap_legend = NULL
  }
  
  root <- rootnode(tree)
  
  # Draw tree with bootstrap nodes, color and shape mappings, and highlighted nodes
  p <- ggtree(x) + geom_star(aes(x, subset = isTip, starshape = shape, fill = species.x), size = 3, show.legend = FALSE) +
    geom_point2(aes(subset = !isTip & node != root, fill = category), shape = 21, size = 2) + theme_tree(legend.position = c(0.8, 0.2)) +
    scale_fill_manual(values = c(color, bootstrap_colors), guide = bootstrap_legend, name = 'Bootstrap Support (BS)',
                      breaks = c('(75,100]', '(50,75]', '(0,50]'), labels = expression("BS >=75", "50 <= BS < 75", "BS < 50"))
  
  if (!is.null(reference1) || !is.null(reference2)) {
    plot <- p + geom_text(aes(label = ifelse(node %in% ref_species$node, ref_species$name, "")), 
                          position = position_nudge(x = 0.08), vjust = 0.5, hjust = 0.9, size = 3.5, family = "Arial", fontface = "bold")
  } else {
    return(p)
  }
  
  if (!is.null(clades) && !is.null(labels)) {
    plot <- plot + geom_cladelab(node = clades, label = labels, align = TRUE, fill = 'white',
                                 offset.text = labeldist, barsize = 0.9, offset.bar = 0.5, fontsize = fontsize)
  }
  
  return(plot)
}
