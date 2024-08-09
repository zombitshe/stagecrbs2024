######################################
# Load required libraries
library(topGO)
library(ggplot2)
library(DT)
library(circlize)

######################################
### Create TERM2GENE data frame ###
genome_path <- "/enadisk/maison/morlon/stage/data/raw/strigamia-acuminata.b2g.reformated.annot"

# Parse lines into list of vectors (due to strct of .b2g.reformed.annot file)
parsed_lines <- lapply(read_lines(genome_path), function(line) {
  split_line <- strsplit(line, "\t")[[1]]
  if (length(split_line) == 2) {
    # If only 2 columns, add an empty string for the third column
    return(c(split_line, ""))
  } else {
    return(split_line)
  }
})

# Convert list of vectors a data frame, with no description
TERM2GENE <- do.call(rbind, parsed_lines) %>%
  as.data.frame(stringsAsFactors = FALSE) %>% 
  dplyr::select(c(1,2)) %>%                             # no description column
  filter(str_detect(V2, "^GO:")) %>%                     # only lines beginning w/ GO terms
  relocate(V2, .before = V1)                    
colnames(TERM2GENE) <- c("GO_ID", "GENE")               # Annotate data frame


###################################### 
### Create cluster large list (all_cluster_sequences) which maps cluster name to list of genes ###

### Read a txt file and return all seq
ls_gene_cluster <- function(path){
  cluster_genes <- read_lines(path)
  
  cluster_genes <- gsub("[\\\\\"]", "", cluster_genes)
  return(cluster_genes)
}

# read directory
cluster_dir <- list.files("/enadisk/maison/morlon/stage/data/raw/clusters/", full.names = TRUE)
#get name of each file in directory
cluster_name <- list.files("/enadisk/maison/morlon/stage/data/raw/clusters/")
cluster_name <- gsub(".txt", "", cluster_name)
#create large list with all sequences in txt per file
all_cluster_sequences <- lapply(cluster_dir, ls_gene_cluster)
#assign gene cluster name to index elements in large list
all_cluster_sequences <- setNames(all_cluster_sequences, cluster_name)


######################################
# Create a gene to GO mapping (reverse of what we have)
gene2GO <- split(TERM2GENE$GO_ID, TERM2GENE$GENE)

# Get all unique genes
all_genes <- unique(TERM2GENE$GENE)


######################################
### Define function to plot every cluster iteratively
plot_go <- function(gene_cluster, cluster_name, algo_topGO = "weight01", subontology, cutoffpval = 0.01){
  # Create geneList as a named factor
  geneList <- factor(as.integer(all_genes %in% gene_cluster))
  names(geneList) <- all_genes
  
  # Create topGO data object
  GOdata <- new("topGOdata",
                ontology = subontology,  
                allGenes = geneList,
                annot = annFUN.gene2GO, # this maps gene to GO terms, here, anottation is provided as gene-to-GO mapping 
                gene2GO = gene2GO)
  
  # Run enrichment test using elimination method (+classic method)
  resultTopGO.elim <- runTest(GOdata, algorithm = algo_topGO, statistic = "Fisher") # here we choose the algorithm and statistical test we like

  # Generate results table
  allRes <- GenTable(GOdata, elimFisher = resultTopGO.elim, orderBy = "elimFisher", topNodes = 50) # here we should be able to choose which term to keep based on pval...
  #filter based on pval
  allRes$elimFisher <- as.numeric(allRes$elimFisher)
  allRes <- subset(allRes, elimFisher < cutoffpval)
  #################
  # Print the results to the console
  #print(allRes)
  
  
  ### VISUALISATION METHODS ###
  #1: Plot subgraph of significant GO terms
  #showSigOfNodes(GOdata, score(resultTopGO.elim), firstSigNodes = 5, useInfo = 'all')
  
  #2: GeneratePDF of graph
  #printGraph(GOdata, resultTopGO.elim, firstSigNodes = 5, fn.prefix = "myResults", useInfo = "all", pdfSW = TRUE)
  
  #3: scatter plot of p-values
  #plot(score(resultTopGO.classic), score(resultTopGO.elim),
  #     xlab = "classic", ylab = "elim", pch = 20,
  #     main = "Scatter plot of p-values")
  
  #4: interactive table of results
  #datatable(allRes)
    #Annotated: The number of genes in the dataset that are annotated with this GO term.
    #Significant: The number of genes among the significant genes that are annotated with this GO term.
    #Expected: The expected number of significant genes annotated with this GO term based on the proportion of the total number of genes. 

  #5: dot plot
  #pdf(file = paste("/enadisk/maison/morlon/stage/results/topgo_go_enrichment_analysis/", cluster_name, ".pdf", sep = ""))
  #print(
    #ggplot(allRes, aes(x = Term, y = -log10(as.numeric(elimFisher)))) +
    #  geom_point() +
    #  coord_flip() +
    #  theme_minimal() +
    #  labs(title = paste("GO Term Enrichment", cluster_name), 
    #       x = "GO Term", 
    #       y = "-log10(p-value)")
  #)
  #dev.off()
  #################
  return(allRes)  #return the test table to be used in the heatmap generation later on
}


# set name of each cluster
names_clusters <- names(all_cluster_sequences)

# apply function to each cluster of the list + create allRes_list object containing all results for heatmap
allRes_list <- list()
allRes_list <- lapply(seq_along(all_cluster_sequences), function(i) {
  # Call plot_go and store its result (allRes) in the list
  plot_go_result <- plot_go(all_cluster_sequences[[i]], 
                            names_clusters[i],
                            algo_topGO = "elim",##################### Choose "weight", "weight01", "elim", else ######################
                            subontology = "CC", ##################### Choose "BP", "MF", "CC" ######################
                            cutoffpval = 0.001) 
  return(plot_go_result) # Return the allRes table to be added to the list
})
allRes_list <- setNames(allRes_list, names_clusters)




#### DRAW HEATMAP OF ALL CLUSTER ####
# Extract unique terms across all clusters
all_terms <- unique(unlist(lapply(allRes_list, function(x) x$Term)))

# Initialize empty data frame for heatmap
heatmap_data <- data.frame(Term = all_terms)

# Loop each cluster to fill in heatmap data
for (cluster_name in names(allRes_list)) {
  cluster_data <- allRes_list[[cluster_name]]
  
  # Initialize a column for this cluster with default p-value (1)
  heatmap_data[[cluster_name]] <- 1
  
  # Loop through each term in this cluster
  for (term in cluster_data$Term) {
    # Find the row in heatmap_data that matches this term
    row_index <- which(heatmap_data$Term == term)
    
    # Update the p-value for this term in this cluster
    heatmap_data[row_index, cluster_name] <- as.numeric(cluster_data$elimFisher[cluster_data$Term == term])
  }
}

# Remove 'Term' column + convert rest to a matrix
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$Term

# Apply -log10 transformation to the p-values + sort cluster in alphabetical order
heatmap_matrix <- -log10(heatmap_matrix)

#create heatmap
  col_fun = colorRamp2(c(0, max(heatmap_matrix)), c("white", "red"))
  lgd = Legend(col_fun = col_fun, title = "-log(10)p-value")
  Heatmap(heatmap_matrix,
          column_title = "Strigamia Acuminata GO enrichement analysis (CC)",
          column_order =colnames(heatmap_matrix),
          col = col_fun,
          rect_gp = gpar(col = "black", lwd = 0.1),
          show_row_dend = FALSE,
          width = unit(16, "cm"), 
          height = unit(20, "cm"),
          show_heatmap_legend = FALSE)
  draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  
  

  
# TODO
    # put everything into a function which would take as input:
      # the path of the b2g.reformed.annot file
      # the path of the directory in which they're the clusters
      # which heatmap(s) to generate (by default, MF + BP + CC)
      # (generate individual result table ?)
