# Library of functions for constructing multi-level heatmaps using the ComplexHeatmap package.

load_packages <- function( tools ) {
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

load_packages( dependencies )


complex_heatmap <- function (counts, reordered_rows = NULL, reordered_cols = NULL, color_palette = NULL, 
                             head_annotation = NULL, row_annotation = TRUE, gaps_row = NULL, gaps_col= NULL, ...) {
  
    counts <- read.delim(counts)
    counts <- counts[,-ncol(counts)]
    
    ## First column as rownames - useful for horiz_annotation
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
    
    ## Generate color_palette
    col2 <-colorRampPalette(color_palette)(400)
    
    counts_m <- as.matrix.data.frame(counts)
    
    ## Left row annotation for CYP clans ##
    species_ha = rowAnnotation("Species" = anno_text(rownames(counts), gp=gpar(fontsize=11,fontface="italic")))
    
    ## add block annotation for CYP clans ##
    if ( !is.null(head_annotation)) {
      top_annotation = HeatmapAnnotation(clans = anno_block(gp = gpar(fill = head_annotation), 
                                                            labels = names(head_annotation),
                                                            width = unit(0.5,"mm"), labels_gp = gpar(col = "black")))
      }
   
    ## add row annotation with gene counts ##
    if ( row_annotation == TRUE) {
      row_ha = rowAnnotation("Total"=anno_barplot(rowSums(counts_m), border = F,
                                                bar_width = 0.8,
                                                gp = gpar(fill = "azure2",fontsize=40),
                                                add_numbers = T, numbers_rot=0, numbers_offset=unit(1,"mm"),
                                                height=unit(6,"mm"), ylim=c(0,25)))
    }
    
    ## draw heatmap with gene counts ##
    p <- ComplexHeatmap::pheatmap(counts_m, cluster_cols = F, cluster_rows = F, scale="none", 
                                  number_color = "black", gaps_row = gaps_row, gaps_col = gaps_col, 
                                  cellwidth = 17, cellheight = 17, color=col2,
                                  border_color = "white", silent = F, show_colnames = T, show_rownames = F,
                                  display_numbers = F, angle_col = c("45"),
                                  fontsize="15",fontsize_row = 17, fontsize_col = 10,legend = T, 
                                  top_annotation=top_annotation,
                                  left_annotation=species_ha,right_annotation=row_ha)
    
    return (p)
}
