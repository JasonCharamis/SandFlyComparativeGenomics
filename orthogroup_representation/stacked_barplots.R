source("~/bin/plotter.R")

dependencies <- c("ggplot2", "ggthemes","dplyr")
sapply (dependencies, function (y) { package_install(y) })

x <- read.delim("orthology_results_for_R.txt", header = T)

order_vector <- c(118:126, 136:144, 163:171, 73:81, 
                  82:90, 91:99, 127:135, 154:162,
                  109:117, 100:108, 145:153)

x <- x[order_vector,]

species_order <- unique(x$Species)
x$Species <- factor(x$Species, levels = unique_order)  # Set custom levels for the factor

stack_order <- rev(unique(x$Type))[c(1,4:7,3,2,8,9)]
x$Type <- factor(x$Type, levels = stack_order)

ggplot(x, aes(x =Species, y = Number_of_Genes, fill = Type)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = rev(c("darkblue","skyblue2" ,"lightblue2","orange1",
                                   "brown3","salmon","yellow2","green4","grey"))) +
  coord_flip() + theme_few() + scale_color_manual(values="black")

  
  
