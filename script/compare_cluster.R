library(clusterProfiler)
library(dplyr)
library(cowplot)

compare_cluster <- function(rds_folder, result_folder){
#### load rds and combine into a list ####
# Set the folder path where the .rds files are located
folder_path <- rds_folder

# Get the list of .rds files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)

# Initialize an empty list to store the enrich results
enrich_results_list <- list()

# Read and name the enrich results for each .rds file
for (file in file_list) {
  # Extract the cluster name from the file name
  cluster_name <- gsub(".*/|\\.rds", "", file)
  
  # Read the .rds file and assign it to the corresponding cluster name
  enrich_results_list[[cluster_name]] <- readRDS(file)
}
##########################################

#### extract data and create compareClusterobject ####
# Function to extract data from each enrichResult
extract_enrich_data <- function(enrich_result, group_name) {
  data <- enrich_result@result
  data$Cluster <- group_name
  return(data)
}
# Apply the function to all enrichResult objects
all_data <- do.call(rbind, mapply(extract_enrich_data, 
                                  enrich_results_list, 
                                  names(enrich_results_list), 
                                  SIMPLIFY = FALSE))
# Ensure Cluster is a factor with levels in the order of the original list
all_data$Cluster <- factor(all_data$Cluster, levels = names(enrich_results_list))

# Create a list of gene IDs for each cluster
gene_clusters <- lapply(enrich_results_list, function(x) unique(unlist(strsplit(x@result$geneID, "/"))))

# Create a compareClusterResult object
comparison_result <- new("compareClusterResult", 
                         compareClusterResult = all_data,
                         geneClusters = gene_clusters)
##########################################

#### run desired analysis ####
### DOTPLOT ###
pdf(file = paste(result_folder, "compare_cluster_dotplot.pdf"), 
    width = 13,
    height = 25)
dotplot(comparison_result)
dev.off()

### CNETPLOT ###
pdf(file = paste(result_folder, "compare_cluster_cnetplot.pdf"), 
    width = 30,
    height = 30)
cnetplot(comparison_result, 
         node_label = "category", 
         cex.params = list(category_node = 5, 
                           category_label = 1.5),
         layout = "kk")
dev.off()

### P-VALUES PLOT ###
#define pbar function
pbar <- function(x) {
  pi=seq(0, 1, length.out=11)
  
  mutate(x, pp = cut(p.adjust, pi)) %>%
    group_by(pp) %>% 
    summarise(cnt = n()) %>% 
    ggplot(aes(pp, cnt)) + geom_col() + 
    theme_minimal() +
    xlab("p value intervals") +
    ylab("Frequency") + 
    ggtitle("p value distribution")
}

# Run pbar function for each object in enrich_results_list
p_list <- lapply(enrich_results_list, pbar)

# save all pbar results
pdf(file = paste(result_folder, "clusters_pvalue_distribution.pdf"), 
    width = 15,
    height = 15)
cowplot::plot_grid(plotlist = p_list, ncol = 2, labels = names(enrich_results_list), label_x = 0.5)
dev.off()
}
##########################################

compare_cluster(rds_folder = "/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/",
result_folder = "/enadisk/maison/morlon/stage/results/compare_cluster/")