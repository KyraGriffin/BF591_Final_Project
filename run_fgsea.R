library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')



#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  
  colnames(deseq_data)[1] <- "genes"
  labeled <- deseq_data %>% 
    as_tibble() %>% 
    mutate(volc_plot_status = case_when(log2FoldChange > 0 & padj < padj_threshold ~ 'UP', 
                                        log2FoldChange < 0 & padj < padj_threshold ~ 'DOWN', 
                                        TRUE ~ 'NS'))
  return(labeled)
}

#' Function to run fgsea on DESeq2 results
#'
#' @param labeled_results (tibble): the labeled results from DESeq2
#' @param gmt (str): the path to the GMT file
#' @param min_size: the threshold for minimum size of the gene set
#' @param max_size: the threshold for maximum size of the gene set
#'
#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
run_gsea <- function(labeled_results, min_size, max_size) {
  
  labeled_results <- labeled_results %>% separate(genes, sep='\\.', into='genes', remove=TRUE)
  gene_ids <- labeled_results %>% pull(genes)
  
  # human <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  # mouse <- useMart('ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
  # 
  # hgnc_symbols <- getLDS(attributes=c('ensembl_gene_id'), 
  #                        filters='ensembl_gene_id', 
  #                        values=gene_ids, 
  #                        mart = mouse, 
  #                        attributesL = c('hgnc_symbol'), 
  #                        martL = human, 
  #                        uniqueRows= TRUE)
  
  hgnc_symbols <- read_csv("data/biomart_ensmusg_to_hgnc.csv", show_col_types = FALSE)
  
  hgnc_results <- labeled_results %>% left_join(hgnc_symbols, by=c('symbol' = 'HGNC.symbol'))
  
  rnks <- hgnc_results %>% 
    drop_na(symbol, log2FoldChange) %>% 
    distinct(symbol, log2FoldChange, .keep_all=TRUE) %>%
    arrange(desc(log2FoldChange)) %>% 
    dplyr::select(symbol, log2FoldChange) %>% 
    deframe()
  
  c2_pathways <- gmtPathways("data/c2.cp.v7.5.1.symbols.gmt")
  
  fgsea_results <- fgsea(c2_pathways, rnks, 
                         minSize = min_size, maxSize = max_size) %>% as_tibble()
  
  return(fgsea_results)
}


########## FUNCTION CALL ##############

deseq_data <- read_tsv("data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv", 
                       col_names = TRUE,
                       show_col_types = FALSE)

labeled_results <- label_res(deseq_data, 0.10) 

labeled_results %>% 
  arrange(padj) %>% 
  relocate(genes, volc_plot_status, log2FoldChange, padj)

########## FUNCTION CALL ##############

fgsea_results <- run_gsea(labeled_results, 15, 500)
fgsea_results
fgsea_results <- write_csv(fgsea_results, "data/fgsea.csv")

####################

####################