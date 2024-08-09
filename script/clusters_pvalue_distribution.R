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

# Plot all pbar results
# Now save cnetplot
pdf(file = "/enadisk/maison/morlon/stage/results/clusters_pvalue_distribution/clusters_pvalue_distribution.pdf", 
    width = 15,
    height = 15)
cowplot::plot_grid(plotlist = p_list, ncol = 2, labels = names(enrich_results_list), label_x = 0.5)
dev.off()
