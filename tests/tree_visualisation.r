library(ggplot2)
library(Biostrings)
library(ggtree)
library("ape")
library(phytools)

# Import tree and csv as param arguments

# Import phylogenetic Tree
tree<- read.tree('/home/phe.gov.uk/nicholas.ellaby/Documents/git_repos/dev_code/usher_genotyping/test_data/usher_out/final-tree.nh')
# Import annotated_summary_stats.csv
meta<-read.csv('/home/phe.gov.uk/nicholas.ellaby/Documents/git_repos/dev_code/usher_genotyping/test_data/annotated_summary_stats.csv', header=TRUE)
# make sure tip.labels match 
row.names(meta)<- meta$sample_id
# Root Tree

# Colour palettes
## Palette for unique genotypes

## Palettes for change in genotypes


# Generate Tree
p <- ggtree(tree) %<+% meta + geom_tippoint(aes(color=original_clade)) # exploring original assignments
p2 <- ggtree(tree) %<+% meta + geom_tippoint(aes(color=estimated_clade)) # exploring new assignments
p3 <- ggtree(tree) %<+% meta + geom_tippoint(aes(color=quality_group)) # annotate quality of tree assignment 


# Subset relevant information for heatplot
heaplot_subset_1 <- meta[,c("original_clade","estimated_clade","clade_assignment_summary")]

# Generate heatplots alongside the tree
gheatmap(p3, heaplot_subset_1, offset = 0, width=0.1, font.size=3, colnames_angle=-45, hjust=0)+theme(legend.position = c(0.075, 0.75))
