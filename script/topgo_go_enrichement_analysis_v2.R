######################################
# Load required libraries
library(topGO)
library(ggplot2)
library(DT)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(stringr)

######################################

go_analysis <- function(b2g, cluster_directory, algo = "weight01", subGO = "BP", cutoffpvalue = 0.05, organism_name = "sample", cutoffclusterterms = 50){
  
  #### Create TERM2GENE data frame ####
  genome_path <- b2g
  
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
  
  ###################################
  
  #### Create cluster large list (all_cluster_sequences) which maps cluster name to list of genes ####
  
  ### Read a txt file and return all seq
  ls_gene_cluster <- function(path){
    cluster_genes <- read_lines(path)
    
    cluster_genes <- gsub("[\\\\\"]", "", cluster_genes)
    return(cluster_genes)
  }
  
  # read directory
  cluster_dir <- list.files(cluster_directory, full.names = TRUE)
  #get name of each file in directory
  cluster_name <- list.files(cluster_directory)
  cluster_name <- gsub(".txt", "", cluster_name)
  #create large list with all sequences in txt per file
  all_cluster_sequences <- lapply(cluster_dir, ls_gene_cluster)
  #assign gene cluster name to index elements in large list
  all_cluster_sequences <- setNames(all_cluster_sequences, cluster_name)
  ###################################
  
  #### GENE2GO + get genes ####
  # Create a gene to GO mapping (reverse of what we have)
  gene2GO <- split(TERM2GENE$GO_ID, TERM2GENE$GENE)
  
  # Get all unique genes
  all_genes <- unique(TERM2GENE$GENE)
  
  # set name of each cluster
  names_clusters <- names(all_cluster_sequences)
  
  ###################################
  
  #### Define function to plot every cluster iteratively #########
  plot_go <- function(gene_cluster, cluster_name, algo_topGO = "weight01", subontology, cutoffpval = 0.01){
    # Create geneList as a named factor
    geneList <- factor(as.integer(all_genes %in% gene_cluster))
    names(geneList) <- all_genes
    
    # Create topGO data object
    GOdata <- new("topGOdata",
                  ontology = subontology,  
                  allGenes = geneList,
                  annot = annFUN.gene2GO, # this maps gene to GO terms, here, annotation is provided as gene-to-GO mapping 
                  gene2GO = gene2GO)
    
    # Run enrichment test using elimination method (+classic method)
    resultTopGO.elim <- runTest(GOdata, algorithm = algo_topGO, statistic = "Fisher") # here we choose the algorithm and statistical test we like
    
    # Generate results table
    allRes <- GenTable(GOdata, elimFisher = resultTopGO.elim, orderBy = "elimFisher", topNodes = cutoffclusterterms) # here we should be able to choose which term to keep based on pval...
    #filter based on pval
    allRes$elimFisher <- as.numeric(allRes$elimFisher)
    allRes <- subset(allRes, elimFisher < cutoffpval)

    return(allRes)  #return the test table to be used in the heatmap generation later on
  }
  ###################################
  
  #### apply function to each cluster of the list + create allRes_list object containing all results for heatmap ####
  allRes_list <- list()
  allRes_list <- lapply(seq_along(all_cluster_sequences), function(i) {
    # Call plot_go and store its result (allRes) in the list
    plot_go_result <- plot_go(all_cluster_sequences[[i]], 
                              names_clusters[i],
                              algo_topGO = algo,##################### Choose "weight", "weight01", "elim", else ######################
                              subontology = subGO, ##################### Choose "BP", "MF", "CC" ######################
                              cutoffpval = cutoffpvalue) 
    return(plot_go_result) # Return the allRes table to be added to the list
  })
  allRes_list <- setNames(allRes_list, names_clusters)
  ###################################
  
  #### setup heatmap for all clusters ####
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
  ###################################

  #### create heatmap ####
  col_fun = colorRamp2(c(0, max(heatmap_matrix)), c("white", "red"))

  lgd = Legend(col_fun = col_fun, title = "-log(10)p-value")

  ht <- Heatmap(heatmap_matrix,
          column_title = paste(organism_name, (" GO enrichement analysis "), subGO, sep = ""),
          column_order =colnames(heatmap_matrix),
          col = col_fun,
          rect_gp = gpar(col = "black", lwd = 0.1),
          show_row_dend = FALSE,
          width = unit(16, "cm"), 
          height = unit(20, "cm"),
          show_heatmap_legend = FALSE)
  draw(ht)
  draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  
  ###################################
}


go_analysis(b2g = "/enadisk/maison/morlon/stage/data/raw/strigamia-acuminata.b2g.reformated.annot", # input b2g.reformated.annot file
            cluster_directory = "/enadisk/maison/morlon/stage/data/raw/clusters/",                  # input clusters directory (contains .txt)
            algo = "weight01",                                                                      # choose topGO algorithm (classic, elim, weight, weight01 (default),...)
            subGO = "BP",                                                                           # choose GO sub-ontology for enrichment analysis
            cutoffpvalue = 0.005,                                                                   # define pvalue cutoff 
            organism_name = "Strigamia Acuminata",                                                   # define organism/sample name for final heatmap
            cutoffclusterterms = 50
            )
