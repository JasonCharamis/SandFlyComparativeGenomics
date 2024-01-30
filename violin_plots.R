
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

dependencies <- c("svglite", "ggplot2", "ggrepel", "ggthemes", "stringr", "dplyr", "ggbreak","tidyverse")

for ( pkg in dependencies ) {
  package_install(pkg)
}

## Function to color points in scatterplot, including distinct colors for specified gene ids
mark_genes <- function(expression_df, selected_genes, selected_genes_color, parent_color, transparency) { 
  
  expression_df$PointColor <- "snow4"
  
  if (is.null(selected_genes)) {
    expression_df <- expression_df %>% mutate(PointColor = "black")
  } else {
    expression_df <- expression_df %>% mutate(PointColor = sapply(GeneID, function(x) { 
      if (grepl(x, selected_genes)) { 
        PointColor = selected_genes_color
      } else {
        PointColor = alpha(parent_color, transparency)
      }
    }))
  }
  return(expression_df)
}


## Function to create violin plots using the mean of gene expression data 
expression_violin_plot <- function (filename, title = "Violin Plot of Gene Expression", selected_genes = NULL, 
                                    mode = "mean", selected_genes_color = "red", 
                                    parent_color = "grey", transparency = 0.5, ... ) {
  
  ## open file and keep Species name ##
  expression_data <- read.delim(filename, header = FALSE)
  species_name <- gsub("_.*", "", expression_data[, 1])
  expression_data$Species <- species_name
  expression_data <- data.frame(GeneID = expression_data[, 1], expression_data[, 2:ncol(expression_data)])
  
  ## Dynamically identify expression column names, while also avoiding the first (GeneID) and last (Mean_Expression-will be added) columns
  columns_to_average <- names(expression_data)[-c(1, ncol(expression_data))]
  
  if ( mode == "mean" ) {
    for (i in 1:nrow(expression_data)) {
      expression_data$Mean_Expression[i] <- rowMeans(expression_data[i, columns_to_average], na.rm = TRUE)
    }
  }
  
  else if ( mode == "median" ) {
    for (i in 1:nrow(expression_data)) {
      expression_data$Median_Expression[i] <- rowMedians(expression_data[i, columns_to_average], na.rm = TRUE)
    }
  }
  
  else {
    print ("Pick a valid option: Mean or Median value of replicates.")
  }
 
  if ("Mean_Expression" %in% colnames(expression_data) || "Median_Expression" %in% colnames(expression_data)) {
          selected_genes_df <- read.delim(selected_genes, header = FALSE)
          expression_data <- mark_genes(expression_df = expression_data, selected_genes = selected_genes_df, 
                                        selected_genes_color = selected_genes_color, parent_color = parent_color, 
                                        transparency = transparency)
  }
  
  if (nrow(expression_data) == 0) {
    stop("Error: No data available after processing.")
  } else { ## Use ggplot2 with geom_violin and geom_jitter to create the violin plot
        ggplot(expression_data, aes(x = Species, y = Mean_Expression)) +
          geom_violin(trim = FALSE, fill = "grey", draw_quantiles = c(0, 0.25, 0.5, 0.75), scale = "area", alpha = 0.05) + 
          geom_jitter(aes(color = PointColor, fill = PointColor), width = 0.4, size = 2) +
          scale_color_identity() + scale_fill_identity() + 
          labs(title = title, y = "Normalized Expression") + 
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
  }
}