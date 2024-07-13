#!/usr/bin/env Rscript

install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, Biostrings, ggtree, ape, phytools, ggnewscale, ggstar, ggtreeExtra)

# # Import tree and csv as param arguments
args = commandArgs(trailingOnly=TRUE)

# # test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("three arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  tree_fp = args[1]
  metadata_fp = args[2]
  output_dir_fp = args[3]
}

# # Import phylogenetic Tree
tree<- read.tree(tree_fp)
# Import annotated_summary_stats.csv
meta<-read.csv(metadata_fp, header=TRUE)
# make sure tip.labels match 
row.names(meta)<- meta$sample_id
# # Root Tree

# Each Genotyping Tree
p1_fn = paste(output_dir_fp, 'tree_with_original_clade_designations.png', sep='/')
# Generate Tree
p <- ggtree(tree) %<+% meta + geom_tippoint(aes(color=original_clade), size=0.5) # exploring original assignments
ggsave(p1_fn, p, height=20, width=15, units='cm', dpi=400)

p2_fn = paste(output_dir_fp, 'tree_with_esimated_clade_designations.png', sep='/')
p2 <- ggtree(tree) %<+% meta + geom_tippoint(aes(color=estimated_clade), size=0.5) # exploring new assignments
ggsave(p2_fn, p2, height=20, width=15, units='cm', dpi=400)

# Placement quality tree
p3_fn = paste(output_dir_fp, 'tree_with_esimated_placement_confidence.png', sep='/')
p3 <- ggtree(tree,) %<+% meta + geom_tippoint(aes(color=quality_group), size=0.5) + scale_colour_discrete(name="Placement Quality") # annotate quality of tree assignment 
ggsave(p3_fn, p3, height=20, width=15, units='cm', dpi=400)


# Subset relevant information for heatplot
heatplot_subset_1 <- meta[,c("original_clade","estimated_clade")]
names(heatplot_subset_1)<-c("Original", "Estimated")

# Summary Tree
p6_fn = paste(output_dir_fp, 'Summary_tree_with_heatplot.png', sep='/')
p4 <- gheatmap(p3, heatplot_subset_1, offset = 0, width=0.1, font.size=3, colnames_angle=-90, hjust=0) + scale_fill_viridis_d(name="Genotypes")

heatplot_subset_2 <- as.data.frame(meta[,"clade_assignment_summary"])
names(heatplot_subset_2) <- c("Assignment Status")

row.names(heatplot_subset_2) <- row.names(meta)
p5 <- p4 + new_scale_fill()

p6 <- gheatmap(p5, heatplot_subset_2, offset=20, width=0.05, font.size=3, colnames_angle=-90, hjust=0 ) + scale_fill_brewer(palette = "Dark2", name="Assignment Status") +
    ggtree::vexpand(.1, -1) # Expand to ensure labels aren't cut-off

ggsave(p6_fn, p6, height=20, width=15, units='cm', dpi=400)
