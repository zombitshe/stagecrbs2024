library(clusterProfiler) 
library(readr)
library(stringr)
library(dplyr)
library(GO.db)
library(ggplot2)
library(cowplot)

#################################################################################################
### Create TERM2GENE data frame
# Define genome genes path
genome_path <- "/enadisk/maison/morlon/stage/data/raw/strigamia-acuminata.b2g.reformated.annot"

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
  filter(str_detect(V1, "^strigamia")) %>%              # only lines w/ strigamia genes
  filter(str_detect(V2, "GO:")) %>%                     # only lines beginning w/ GO terms
  relocate(V2, .before = V1)                    
colnames(TERM2GENE) <- c("GO_ID", "GENE")               # Annotate data frame

#################################################################################################
### TERM2NAME, link between GO_ID and GO_TERM using GO.db package
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

#################################################################################################
### Create gene cluster list
# Define path
cluster_path <- "/enadisk/maison/morlon/stage_martin/stage24/results/gene_clusterH.txt"       # CHANGE SOURCE CLUSTER HERE - STAGE MARTIN
# Create list 
cluster_genes <- grep("strigamia-acuminata", read_lines(cluster_path), value = TRUE)

cluster_genes <- gsub("[\\\\\"]", "", cluster_genes)  #delete residual " and / characters

#################################################################################################
### Perform enrichment analysis & save it
enrich_result <- enricher(gene = cluster_genes,
                          TERM2GENE = TERM2GENE,
                          TERM2NAME = TERM2NAME,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 1,
                          maxGSSize = 1000) 
#saveRDS(enrich_result, file = "/enadisk/maison/morlon/stage/results/go_enrichment_analysis/GO_analysis_strigamia_acuminata_cluster_H.rds")

# Draw plot & save it
#svg("/enadisk/maison/morlon/stage/results/go_enrichment_analysis/GO_analysis_strigamia_acuminata_cluster_H.svg")
barplot(enrich_result) + ggtitle("GO enrichment analysis - Cluster H Strigamia Acuminata", 
                                 subtitle = paste("Cluster size =", length(cluster_genes)))
#dev.off()



##########################
### TODO
#~~Add cluster size~~
# separate in two functions, script should read a folder of GO clusters
# explore the possibility to use the enrichGO() function instead, how does this affect the results ? 
  # with GO.db as OrgDb annotation database argument (https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData)
  # a notable advantage would be the control of which type of GO term to search for during enrichment analysis + control of pooling of GO sub-ontology
