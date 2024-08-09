#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#install.packages("amap")
#install.packages('dynamicTreeCut')
#install.packages('phylogram')
#install.packages("dendextend")

library(ComplexHeatmap)         # generate heatmaps
library(amap)                   # compute distances
library(dynamicTreeCut)         # dendrogram cut
library(stats)
library(ggtree)                 # plot phylogenetic tree
library(ape)                    # manipulate phylogenetic tree
library(dplyr)
library(phylogram)
library(dendextend)             # manipulate dendrogram
library(ggplot2)
library(colorspace)

### Paths definition
phylo_file <- "/gstock/MetaInvert/martin/tree_ortho_eukaryota_myriapoda.nwk"
profiles_file <- '/gstock/MetaInvert/martin/strigamia-acuminata_profiles.tsv.short'
info_species_file <- "/gstock/MetaInvert/martin/info_about_species.csv"

### Phylo tree processing
# Read the Newick formatted tree
phylo_tree <- read.tree(phylo_file)
dendrogram_phylo <- as.dendrogram.phylo(phylo_tree)

# save dendrogram
#svg("/home/merlat/projects/myriapods/orthoinspector/results/phylo_dendrogram.svg", width = 30 ,height = 32)
#par(cex=0.1, mar=c(5, 8, 4, 1))
#plot(dendrogram_phylo, xlab="", ylab="", main="", sub="", axes=FALSE)
#dev.off()

### Profiles analysis
# Load the the profiles
profiles_matrix <- as.matrix(read.table(profiles_file, header = TRUE, row.names = 1, sep = "\t", check.names=FALSE))

# Calculate distances between profiles
distance <- Dist(profiles_matrix, method = "euclidean")

# Hierarchical clustering
my_hclust_gene <- hclust(distance, method = "ward.D2")

#~ # Cut the tree to get clusters
#~ my_gene_col <- cutreeDynamic(dendro = my_hclust_gene, method = "hybrid", deepSplit = 4, distM = as.matrix(distance))                                                       # A modifier

#~ # Merge cluster information with profiles
#~ df <- cbind(my_gene_col, profiles_matrix)
#~ df <- df[order(my_gene_col),] 
#~ matrix <- data.matrix(df[,-1])

# Generate dendrogram
dendrogram_clustering <- as.dendrogram.phylo(my_hclust_gene)

### Annotation
# Define the factor_id function
factor_id <- function(df, phylo) {
  # Set taxid as a factor with levels from phylo$tip.label
  levels <- phylo$tip.label
  df$taxid <- factor(df$taxid, levels = levels)
  
  # Sort the dataframe by taxid
  df <- df %>% arrange(taxid)
  
  return(df)
}

taxo <-  read.csv(file = info_species_file)
taxo <- factor_id(taxo, dendrogram_phylo)

# # Define the list of clades of interest
clades_of_interest <- c("Arthropoda")
                        #,"Ascomycota", "Chordata", "Streptophyta", 
                        #"Platyhelminthes", "Microsporidia", "Oomycota", "Mucoromycota", 
                        #"Placozoa", "Apicomplexa", "Basidiomycota", "Nematoda", 
                        #"Rotifera", "Evosea", "Chytridiomycota", "Parabasalia", 
                        #"Chlorophyta", "Ciliophora", "Cryptomycota", "Discosea", 
                        #"Rhodophyta", "Zoopagomycota", "Euglenozoa", "Echinodermata", 
                        #"Haptophyta", "Bacillariophyta", "Fornicata", "Cnidaria", 
                        #"Annelida", "Mollusca", "Perkinsozoa", "Endomyxa", "Porifera", 
                        #"Foraminifera", "Heterolobosea", "Blastocladiomycota", "Brachiopoda")
##### Define Subphylium + Class of interest
subphylium_of_interest <- c("Myriapoda")
class_of_interest <- c("Chilopoda", "Diplopoda")
# 
# # Replace clades not in the list with NA
taxo <- taxo %>%
  mutate(phylum_taxid = ifelse(phylum_taxid %in% clades_of_interest, phylum_taxid, NA))
#### Same with subphylium/class
taxo <- taxo %>%
  mutate(subphylum_taxid = ifelse(subphylum_taxid %in% subphylium_of_interest, subphylum_taxid, NA))

taxo <- taxo %>%
  mutate(class_taxid = ifelse(class_taxid %in% class_of_interest, class_taxid, NA))

# define color palettes for clades with unique colors
colors <- colorspace::rainbow_hcl(length(clades_of_interest))
clade_colors <- setNames(colors, clades_of_interest)
#### define color palettes for subphylium/class with unique colors
colors <- colorspace::rainbow_hcl(length(subphylium_of_interest))
subphylium_colors <- setNames(colors, subphylium_of_interest)

colors <- colorspace::rainbow_hcl(length(class_of_interest))
class_colors <- setNames(colors, class_of_interest)

#### Arrange data by groups
taxo <- arrange(taxo, kingdom_taxid, phylum_taxid, subphylum_taxid, class_taxid)
#### Annotate heatmap
taxo_annotation <- HeatmapAnnotation(na_col = "white", 
                                     Kingdom = taxo$kingdom_taxid, 
                                     Phylum = taxo$phylum_taxid, 
                                     Subphylium = taxo$subphylum_taxid,
                                     Class = taxo$class_taxid,
                                     col = list(Kingdom = c(Metazoa = "red", Fungi = "blue", Viridiplantae = "green"), 
                                                Phylum = clade_colors, 
                                                Subphylium = subphylium_colors,
                                                Class = class_colors))
#install.packages("dendsort")
#library(dendsort)
#row_dend = dendsort(hclust(dist(profiles_matrix)))
#col_dend = dendsort(hclust(dist(t(profiles_matrix))))
# Generate the heatmap
ht <- Heatmap(profiles_matrix, cluster_columns = dendrogram_phylo, clustering_method_row = "ward.D2", 
              clustering_distance_rows = "euclidean",
              show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = TRUE,
              show_column_dend = FALSE,
              row_dend_width = unit(2, "cm"), column_dend_height = unit(2, 'cm'),
              top_annotation = taxo_annotation)
ht
# Download the heatmap
#svg("/home/schoenstein/stage24/heatAnno/heatmapEuclideanComplete.svg", width = 30 ,height = 32)
#ht = Heatmap(profiles_matrix, cluster_columns = dendrogram_phylo, clustering_method_row = "ward.D2", 
#        clustering_distance_rows = "euclidean",
#        show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = TRUE,
#        row_dend_width = unit(2, "cm"), column_dend_height = unit(2, 'cm'),
#        top_annotation = taxo_annotation)
#draw(ht)
#dev.off()




