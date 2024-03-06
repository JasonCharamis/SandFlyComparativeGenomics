library(BiocManager)
library(ggplot2)

package_install <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    print(sprintf("%s %s", package_name, "is not installed. Installing it!"))
    
    if (BiocManager::available(package_name)) {
      BiocManager::install(package_name, dependencies = TRUE, update = TRUE)
    } else {
      install.packages(package_name, dependencies = TRUE, ask = FALSE, reinstall = TRUE)
    }
  }
  library(package_name, character.only = TRUE)
}

dependencies <- c("svglite", "ggplot2", "ggrepel", "ggthemes", "stringr", "dplyr", "ggpubr", "tidyverse")

for (pkg in dependencies) {
  package_install(pkg)
}

.mark_genes <- function(expression_df, selected_genes, parent_color = "grey") {
  expression_data$PointColor <- parent_color
  
  for (id in names(selected_genes)) {
    i <- which(expression_data$GeneID == id)
    expression_data$PointColor[i] <- selected_genes[id]
  } 
}

expression_violin_plot <- function(filename, title = "Violin Plot of Gene Expression", selected_genes = NULL,
                                   mode = "mean", parent_color = "grey", transparency = 0.5, ...) {
  
  expression_data <- read.delim(filename, header = FALSE)
  species_name <- gsub("_.*", "", expression_data[, 1])
  expression_data$Species <- species_name
  expression_data <- data.frame(GeneID = expression_data[, 1],
                                expression_data[, 2:ncol(expression_data)])
  
  columns_to_average <- names(expression_data)[-c(1, ncol(expression_data))]
  
  if (mode == "mean") {
    expression_data$Mean_Expression <- rowMeans(expression_data[, columns_to_average], na.rm = TRUE)
  } else if (mode == "median") {
    expression_data$Median_Expression <- apply(expression_data[, columns_to_average], 1, median, na.rm = TRUE)
  } else {
    stop("Pick a valid option: Mean or Median value of replicates.")
  }
  
  if (!is.null(selected_genes)) {
      expression_data <- .mark_genes(expression_df = expression_data, 
                                     selected_genes = selected_genes,
                                     parent_color = parent_color)
  }
  
  ggplot(expression_data, aes(x = Species, y = Mean_Expression)) +
  geom_violin(trim = FALSE, fill = parent_color, draw_quantiles = c(0, 0.25, 0.5, 0.75), scale = "area", alpha = 0.05) +
  geom_jitter(aes(fill = PointColor), width = 0.4, size = 2) +
  scale_color_identity() + scale_fill_identity() +
  labs(title = title, y = "Normalized Expression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_log10()
}

