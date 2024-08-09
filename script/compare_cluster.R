library(clusterProfiler)
library(dplyr)
library(cowplot)

#load each cluster results
clusterA <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterA.rds")
clusterB <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterB.rds")
clusterC <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterC.rds")
clusterC1 <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterC1.rds")
clusterC2 <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterC2.rds")
clusterD <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterD.rds")
clusterE <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterE.rds")
clusterF <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterF.rds")
clusterG <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterG.rds")
clusterH <- readRDS("/enadisk/maison/morlon/stage/results/go_enrichment_analysis_v2/GO_analysis_s_acuminata_gene_clusterH.rds")

#combine results into list
enrich_results_list <- list(clusterA,
                            clusterB, 
                            clusterC, 
                            clusterC1,
                            clusterC2,
                            clusterD,
                            clusterE,
                            clusterF,
                            clusterG,
                            clusterH)
enrich_results_list_names <- c('clusterA','clusterB','clusterC','clusterC1','clusterC2','clusterD','clusterE','clusterF','clusterG','clusterH')

enrich_results_list <- setNames(enrich_results_list, enrich_results_list_names)


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

# Now save dotplot
pdf(file = "/enadisk/maison/morlon/stage/results/compare_cluster/compare_cluster_dotplot.pdf", 
    width = 13,
    height = 25)
dotplot(comparison_result)
dev.off()

# Now save cnetplot
pdf(file = "/enadisk/maison/morlon/stage/results/compare_cluster/compare_cluster_cnetplot.pdf", 
    width = 30,
    height = 30)
cnetplot(comparison_result, 
         node_label = "category", 
         cex.params = list(category_node = 5, 
                           category_label = 1.5),
         layout = "kk")
dev.off()

### P-VALUES PLOT ##
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
pdf(file = "/enadisk/maison/morlon/stage/results/compare_cluster/clusters_pvalue_distribution.pdf", 
    width = 15,
    height = 15)
cowplot::plot_grid(plotlist = p_list, ncol = 2, labels = names(enrich_results_list), label_x = 0.5)
dev.off()
