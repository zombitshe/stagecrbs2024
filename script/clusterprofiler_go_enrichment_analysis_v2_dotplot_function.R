### INSTRUCTIONS
# input you background gene .b2g.reformed.annot           --> line 14
# input your repertory with all cluster.txt               --> line 66 and 68
# input the repertory where the .rds and .svg are saved   --> line 92 and 95
# run the script

library(clusterProfiler) 
library(readr)
library(stringr)
library(dplyr)
library(GO.db)
library(ggplot2)
library(cowplot)

profiler_go_analysis <- function(b2g, cluster_directory, output_directory, 
plot_analysis = TRUE, cutoffpvalue = 0.05, adjustpvalue= "BH", orga_name = "_s_acuminata_"){
  #### Create TERM2GENE data frame ####
  # Define genome genes path
  genome_path <- b2g

  # Parse lines into list of vectors
  parsed_lines <- lapply(read_lines(genome_path), function(line) {
    split_line <- strsplit(line, "\t")[[1]]
    if (length(split_line) == 2) {
      # If only 2 columns, add an empty string for the third column
      return(c(split_line, ""))
    } else {
      return(split_line)
    }
  })

  # Convert the list of vectors to a data frame, with no description
  TERM2GENE <- do.call(rbind, parsed_lines) %>%
    as.data.frame(stringsAsFactors = FALSE) %>% 
    dplyr::select(c(1,2)) %>%                             # no description column
    #filter(str_detect(V1, "^strigamia")) %>%              
    filter(str_detect(V2, "^GO:")) %>%                     # only lines beginning w/ GO terms
    relocate(V2, .before = V1)                    
  colnames(TERM2GENE) <- c("GO_ID", "GENE")               # Annotate data frame

  ####################################

  #### TERM2NAME, link between GO_ID and GO_TERM using GO.db package ####
  # Create function to associate GO_TERM to GO_ID in a TERM2GENE table, careful to return NA if no TERM found
  get_go_term_safe <- function(go_id) {
    go_term <- tryCatch({
      term <- GOTERM[[go_id]]
      if (is.null(term)) return(NA)  # Return NA if GO ID not found
      Term(term)
    }, error = function(e) {
      NA  # Return NA on error
    })
    return(go_term)
  }

  # create TERM2NAME dataframe
  TERM2NAME <- transform(TERM2GENE, GO_TERM = sapply(GO_ID, get_go_term_safe)) %>% 
    dplyr::select(c(1,3))

  ####################################

  #### Read cluster directory to create large list of all clusters ####
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
  ####################################

  #### execute analysis ####
  # set organism name (for output file)
  org_name <- orga_name

  #for all sequence in the input directory
  for (i in seq_along(all_cluster_sequences)){
    clusters <- all_cluster_sequences[[i]]
    # perform enrichment
    enrich_result <- enricher(gene = clusters,
                              TERM2GENE = TERM2GENE,
                              TERM2NAME = TERM2NAME,
                              pvalueCutoff = cutoffpvalue,
                              pAdjustMethod = adjustpvalue,
                              minGSSize = 1,
                              maxGSSize = 1000)
    # save result as .rds
    saveRDS(enrich_result, 
            paste(output_directory,"GO_analysis", org_name, names(all_cluster_sequences)[i], ".rds", sep = ""))
    # save bar plot
    if (plot_analysis){    
    tryCatch({
      svg_file <- paste(output_directory,"GO_analysis", org_name, names(all_cluster_sequences)[i], ".svg", sep = "")
      svg(svg_file)
      print(dotplot(enrich_result) + ggtitle(paste("enrichGO", org_name, names(all_cluster_sequences)[i]), 
                                            subtitle = paste("Cluster size =", length(clusters))))
      dev.off()
    }, error = function(e) {
      print(paste("Error occurred while drawing SVG for", names(all_cluster_sequences)[i]))
      file.remove(svg_file)  # Remove the SVG file if there is an error
    })
  }
}
}
####################################

profiler_go_analysis(b2g = "/enadisk/maison/morlon/stage/data/raw/strigamia-acuminata.b2g.reformated.annot", 
cluster_directory = "/enadisk/maison/morlon/stage/data/raw/clusters/",
output_directory = "/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2_dotplot/",
plot_analysis = TRUE,
cutoffpvalue = 0.05,
adjustpvalue= "BH",
orga_name = "_s_acuminata_")
