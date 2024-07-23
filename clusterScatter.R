#cluster_result: the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#dimred: dimension reduction method. Options are "UMAP", "tSNE", or "PCoA".
#scores: How should points be sized on the scatterplot? Options are "pvalue" and "size". Pvalue will scale points according to -log10(pvalue), size will scale points according to gene set size K


clusterScatter <- function(
    cluster_result = NULL,
    dimred = "UMAP",
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
  
  olmd <- cluster_result$dist_mat
  
  ## run dimension reduction on 2 axes ##
  if(dimred == "UMAP"){
    custom.settings = umap::umap.defaults
    custom.settings$input = "dist"
    custom.settings$random_state = 432143
    plot_df <- umap::umap(olmd, custom.settings)
    plot_df <- plot_df$layout %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var = "pathway") %>% 
      dplyr::left_join(df, by = c("pathway" = "pathway"))
    
    xlab <- "UMAP1"
    ylab <- "UMAP2"
  
  } else if(dimred == "PCoA"){
    plot_df <- as.data.frame(stats::cmdscale(olmd)) %>% 
      tibble::rownames_to_column(var = "pathway") %>% 
      dplyr::left_join(df, by = c("pathway" = "pathway"))
    
    xlab <- "PCoA1"
    ylab <- "PCoA2"
    
    
  } else if(dimred =="tSNE"){
    plot_df <- Rtsne::Rtsne(olmd, is_distance = TRUE, dims = 2, verbose = F, maxiter = 100)
    plot_df <- as.data.frame(plot_df$Y) %>% 
      dplyr::mutate("pathway" = rownames(olmd), .before = "V1")%>% 
      dplyr::left_join(df, by = c("pathway" = "pathway"))
    
    xlab <- "tSNE1"
    ylab <- "tSNE2"
    
  } else(stop("Valid options for dimred are 'UMAP', 'PCoA', or 'tSNE'."))
  
  ## make plot ##
  plot_df<- plot_df %>% 
    dplyr::mutate(score_scale = scales::rescale(score, to = c(1, 30)),
                  labs = stringr::str_replace_all(pathway, "_", " "),
                  labs = tolower(labs),
                  labs = stringr::str_wrap(labs, width = 25)) 
  
  p <- ggplot2::ggplot(plot_df) +
    ggrepel::geom_text_repel(ggplot2::aes(x = V1, y = V2,label = labs, color = cluster), 
                             max.overlaps = Inf, size = 2, show.legend = F) +
    ggplot2::geom_point(ggplot2::aes(x = V1, y = V2, size = score, 
                                     fill = cluster, color = cluster), 
                        alpha = 0.5, show.legend = T) + 
    ggplot2::geom_text(ggplot2::aes(x = V1, y = V2,color = cluster, label = cluster), 
                       size = 2, show.legend = F) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank()) +
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(ylab) + 
    ggplot2::labs(size = scr_str, color = "Cluster", fill = "Cluster")+
    scale_size_continuous(range = c(1, 10))
  
  return(p)
  
}