#' Create treemap plot(s) of clustered enrichment result
#'
#' @param cluster_result the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#' @param scores How should points be sized on the scatterplot? Options are "pvalue" and "size". Pvalue will scale points according to -log10(pvalue), size will scale points according to gene set size K
#'
#' @return a ggplot object(for one-part enrichments like hypergeometric) or a list of ggplot objects (for two-part like GSEA)
#' @export
#'
#' @examples
#' # Not run
#' # Create res object as described in clusterSets example
#' # clusterTreemap(res)

clusterTreemap <- function(
    cluster_result = NULL,
    scores = "pvalue"
){
  
  . <- cluster <- pathway <- pathway_format <- score <- NULL
  set.seed(432143)
  
  # hypergeometric - one plot
  if(is.null(cluster_result$cluster_membership$sign)){
    df <- cluster_result$input_df
    
    ## format inputs ##
    if(scores == "pvalue"){
      df$score <- -log10(df$pvalue)
      scr_str <- "-log10(p-value)"
    } else if(scores == "size"){
      df$score <- df$K
      scr_str <- "Gene Set Size"
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
  else{

    figlist <- list()
    for(s in unique(cluster_result$cluster_membership$sign)){
      df <- cluster_result$input_df
      
      ## format inputs ##
      if(scores == "pvalue"){
        df$score <- -log10(df$pvalue)
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
