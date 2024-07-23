#NOTE: a lot of this was taken directly from vissE package (https://www.bioconductor.org/packages/release/bioc/html/vissE.html)

#cluster_result: the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#scores: How should points be sized on the scatterplot? Options are "pvalue" and "size". Pvalue will scale points according to -log10(pvalue), size will scale points according to gene set size K

clusterTreemap <- function(
    cluster_result = NULL,
    scores = "pvalue"
){
  
  set.seed(432143)
  
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
    dplyr::mutate(
      pathway_format = tolower(pathway),
      pathway_format = str_replace_all(pathway_format, "_", " "),
      pathway_format = str_replace_all(pathway_format, "-", " "),
      pathway_format = str_replace_all(pathway_format, "/", " "),
      pathway_format = str_replace_all(pathway_format, ":", " "),
      pathway_format = str_replace_all(pathway_format, ",", ""),
      pathway_format = str_replace_all(pathway_format, fixed("["), " "),
      pathway_format = str_replace_all(pathway_format, fixed("]"), " ")) %>%
    dplyr::filter(!is.na(cluster))
  
  p <- ggplot2::ggplot(df3, ggplot2::aes(area = score, fill = as.factor(cluster),
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
