# Library of functions for constructing multi-level heatmaps using the ComplexHeatmap package.

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

# Load required packages or install them if necessary
dependencies <- c("ComplexHeatmap",
            		  "ggplot2",
            		  "dplyr",
            		  "tidyverse")

.load_packages( dependencies )


#' complex_heatmap
#' Function to visualize a heatmap with a variety of metadata and customizable features
#'
#' @param counts Path to the counts file.
#' @param reordered_rows Indices specifying the order of rows. Default is NULL.
#' @param reordered_cols Indices specifying the order of columns. Default is NULL.
#' @param color_palette The color palette for heatmap. Default is NULL.
#' @param top_annotation Annotation for the heatmap header. Default is NULL.
#' @param row_annotation Logical indicating whether to include row annotations. Default is TRUE.
#' @param gaps_row The height of gaps between rows. Default is NULL.
#' @param gaps_col The width of gaps between columns. Default is NULL.
#'
#' @return A heatmap plot constructed using ComplexHeatmap package.
#'
#' @export


complex_heatmap <- function (counts, 
                             reordered_rows = NULL, 
                             reordered_cols = NULL, 
                             color_palette = NULL, 
                             top_annotation = NULL,
                             top_annotation_title = NULL,
                             right_annotation = NULL, 
                             right_annotation_title = "Total",
                             left_annotation = NULL, 
                             left_annotation_title = "Category",
                             gaps_row = NULL, 
                             gaps_col= NULL, ...) {
  
    counts <- read.delim(counts)
    counts <- counts[,-ncol(counts)]
    
    # First column as rownames - useful for horizontal_annotation
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    
    if (!is.null(reordered_rows)) {
      counts <- counts[reordered_rows,]
    }
    
    if (!is.null(reordered_cols)) {
      counts <- counts[,reordered_cols]
    }
  
    counts <- counts[,-2]
    counts <- t(counts)
    counts <- as.data.frame(counts)
    
    # Generate color_palette
    col2 <-colorRampPalette(color_palette)(400)
    
    counts_m <- as.matrix.data.frame(counts)
    
    # Left row annotation 
    if (!is.null(left_annotation)) {
      if (!is.null(top_annotation_title)) {
        left_annotation = rowAnnotation(left_annotation_title = anno_text(rownames(counts),
                                        gp=gpar(fontsize = 1,
                                                 fontface = "italic")
                                       )
        )
      } else {
          stop("Please provide a left_annotation_title.")
      }
    }
    
    # Add block annotation 
    if ( !is.null(top_annotation)) {
      if (!is.null(top_annotation_title)) {
      top_annotation = HeatmapAnnotation(top_annotation_title = anno_block(gp = gpar(fill = top_annotation), 
                                                                           labels = names(top_annotation),
                                                                           width = unit(0.5,"mm"), labels_gp = gpar(col = "black")
                                                                          )
                                         )
    } else {
        stop("Please provide a top_annotation_title.")
    }
  }

    # Add row annotation with gene counts
    if ( !is.null(right_annotation) ) {
      if (!is.null(right_annotation_title)) {
        right_annotation = rowAnnotation(left_annotation_title = anno_barplot(rowSums(counts_m),
                                                                              border = F,
                                                                              bar_width = 0.8,
                                                                              gp = gpar(fill = "azure2",fontsize=40),
                                                                              add_numbers = T, numbers_rot=0, numbers_offset=unit(1,"mm"),
                                                                              height=unit(6,"mm"), ylim=c(0,25)
                                                                             )
                                         )
    } else {
        stop("Please provide a right_annotation_title.")
    }
  }
    
    
    # Draw heatmap with gene counts
    plot <- ComplexHeatmap::pheatmap(counts_m, 
                                  cluster_cols = F,
                                  cluster_rows = F,
                                  scale="none", 
                                  number_color = "black",
                                  gaps_row = gaps_row,
                                  gaps_col = gaps_col, 
                                  cellwidth = 17,
                                  cellheight = 17,
                                  color=col2,
                                  border_color = "white",
                                  silent = F,
                                  show_colnames = T,
                                  show_rownames = F,
                                  display_numbers = F,
                                  angle_col = c("45"),
                                  fontsize="15",
                                  fontsize_row = 17,
                                  fontsize_col = 10,
                                  legend = T, 
                                  top_annotation=top_annotation,
                                  left_annotation=species_ha,
                                  right_annotation=right_annotation
                                  )
    
    if ( exists ("plot")) {
      return (p)
    } else {
      print ("Error. Plot was not generated.")
    }
}
