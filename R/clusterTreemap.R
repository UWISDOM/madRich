#' Create treemap plot(s) of clustered enrichment result
#'
#' @param cluster_result the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#' @param scores How should sections be sized on the treemap? Options are "pvalue" and "size". Pvalue will scale points according to -log10(pvalue), size will scale points according to gene set size K
#'
#' @return a ggplot object(for one-part enrichments like hypergeometric) or a list of ggplot objects (for two-part like GSEA)
#' @export
#'
#' @examples
#' #' # Run enrichment using SEARchways
#' gene_list2 <- list(HRV1 = names(SEARchways::example.gene.list[[1]]),
#'                   HRV2 = names(SEARchways::example.gene.list[[2]]))
#' df1 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'                              category="C5", subcategory="GO:MF", ID="ENSEMBL")
#' df2 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'                               category="C5", subcategory="GO:BP", ID="ENSEMBL")
#' df <- dplyr::bind_rows(df1, df2)
#' res <- clusterSets(df = df, enrich_method="hypergeometric",
#'                    ID = "ENSEMBL",
#'                    collections = c("C5"),
#'                    subcollections = c("C5" = "GO:MF", "C5" = "GO:BP"),
#'                    hclust_height = c(0.7),
#'                    group_name = "HRV1",
#'                    fdr_cutoff = 0.4)
#' clusterTreemap(res)

clusterTreemap <- function(
    cluster_result = NULL,
    scores = "pvalue"
){
  
  . <- cluster <- pathway <- pathway_format <- score <- NULL
  set.seed(432143)
  
  # hypergeometric - one plot
  if(is.null(cluster_result$cluster_membership$sign)){ # if there is no sign column, result is from hypergeometric there is one cluster input
    df <- cluster_result$input_df
    
    ## format inputs ##
    if(scores == "pvalue"){
      df$score <- -log10(df$pval)
      scr_str <- "-log10(p-value)"
    } else if(scores == "size"){
      if(!is.null(df$K)){
        df$score <- df$K
        scr_str <- "Gene Set Size"
      } else if(!is.null(df$size_pathway)){
        df$score <- df$size_pathway
        scr_str <- "Gene Set Size"
      } else{
        df$score <- 1
        print("No valid gene set size columns provided. Please make sure input data has size_pathway or K, or else provide p-values to size points")
        scr_str <- "Gene Set"
      }
    }
    
    
    
    cl_df <- cluster_result$cluster_membership
    df <- df %>% 
      dplyr::left_join(cl_df, by = c("pathway" = "pathway"))
    df <- df %>%
      dplyr::filter(!is.na(cluster)) %>%
      dplyr::mutate(
        pathway_format = tolower(pathway),
        pathway_format = stringr::str_replace_all(pathway_format, "_", " "),
        pathway_format = stringr::str_replace_all(pathway_format, "-", " "),
        pathway_format = stringr::str_replace_all(pathway_format, "/", " "),
        pathway_format = stringr::str_replace_all(pathway_format, ":", " "),
        pathway_format = stringr::str_replace_all(pathway_format, ",", ""),
        pathway_format = stringr::str_replace_all(pathway_format, stringr::fixed("["), " "),
        pathway_format = stringr::str_replace_all(pathway_format, stringr::fixed("]"), " ")) 
    
    
    p <- ggplot2::ggplot(df, ggplot2::aes(area = score, fill = as.factor(cluster),
                                          label = pathway_format, 
                                          subgroup = as.factor(cluster))) +
      treemapify::geom_treemap() +
      treemapify::geom_treemap_subgroup_border(colour = "white", size = 5) +
      treemapify::geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                                             alpha = 0.25, colour = "black",
                                             fontface = "italic") +
      treemapify::geom_treemap_text(colour = "white", place = "centre",
                                    size = 15, grow = TRUE, reflow = TRUE) +
      ggplot2::theme(legend.position = "none")
    
    return(p)
  } 
  # gsea - list of 2 plots
  else{ # if there is a sign column, result is from GSEA and clustering is sign-separated

    figlist <- list()
    for(s in unique(cluster_result$cluster_membership$sign)){
      df <- cluster_result$input_df
      
      ## format inputs ##
      if(scores == "pvalue"){
        df$score <- -log10(df$pval)
        scr_str <- "-log10(p-value)"
      } else if(scores == "size"){
        df$score <- df$K
        scr_str <- "Gene Set Size"
      }
      
      cl_df <- cluster_result$cluster_membership %>% 
        dplyr::filter(sign == s)
      df <- df %>% 
        dplyr::filter(pathway %in% cl_df$pathway) %>% 
        dplyr::left_join(cl_df, by = c("pathway" = "pathway"))%>%
        
        dplyr::filter(!is.na(cluster)) %>%
        dplyr::mutate(
          pathway_format = tolower(pathway),
          pathway_format = stringr::str_replace_all(pathway_format, "_", " "),
          pathway_format = stringr::str_replace_all(pathway_format, "-", " "),
          pathway_format = stringr::str_replace_all(pathway_format, "/", " "),
          pathway_format = stringr::str_replace_all(pathway_format, ":", " "),
          pathway_format = stringr::str_replace_all(pathway_format, ",", ""),
          pathway_format = stringr::str_replace_all(pathway_format, stringr::fixed("["), " "),
          pathway_format = stringr::str_replace_all(pathway_format, stringr::fixed("]"), " ")) 
      
      figlist[[s]] <- local({
        
        . <- pathway_format <- cluster <- score <- pathway <- NULL
          
        df <- df
        s <- s
        ggplot2::ggplot(df, ggplot2::aes(area = score, fill = as.factor(cluster),
                                         label = pathway_format, 
                                         subgroup = as.factor(cluster))) +
          treemapify::geom_treemap() +
          treemapify::geom_treemap_subgroup_border(colour = "white", size = 5) +
          treemapify::geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                                                 alpha = 0.25, colour = "black",
                                                 fontface = "italic") +
          treemapify::geom_treemap_text(colour = "white", place = "centre",
                                        size = 15, grow = TRUE, reflow = TRUE) +
          ggplot2::theme(legend.position = "none")
      })
      
      
    }
    return(figlist)
  }
}
