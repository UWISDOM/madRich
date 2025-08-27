#' Create scatterplot of clustered gene set result
#'
#' @param cluster_result the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#' @param dimred dimension reduction method. Options are "UMAP", "tSNE", or "PCoA".
#' @param scores How should points be sized on the scatterplot? Options are "pvalue" and "size". Pvalue will scale points according to -log10(pvalue), size will scale points according to gene set size K
#'
#' @return a ggplot object(for one-part enrichments like hypergeometric) or a list of ggplot objects (for two-part like GSEA)
#' @export
#'
#' @examples
#' # Run enrichment using SEARchways
#' gene_list2 <- list(HRV1 = names(SEARchways::example.gene.list[[1]]),
#'                   HRV2 = names(SEARchways::example.gene.list[[2]]))
#' df1 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'                              collection="C5", subcollection="GO:MF", ID="ENSEMBL")
#' df2 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'                               collection="C5", subcollection="GO:BP", ID="ENSEMBL")
#' df <- dplyr::bind_rows(df1, df2)
#' res <- clusterSets(df = df, enrich_method="hypergeometric",
#'                    ID = "ENSEMBL",
#'                    collections = c("C5"),
#'                    subcollections = c("C5" = "GO:MF", "C5" = "GO:BP"),
#'                    hclust_height = c(0.7),
#'                    group_name = "HRV1",
#'                    fdr_cutoff = 0.4)
#' clusterScatter(res, dimred = "PCoA", scores = "size")
#' 
clusterScatter <- function(
    cluster_result = NULL,
    dimred = "UMAP",
    scores = "pvalue"
){
  
  . <- plot_df <- df <- scr_str <- cl_df <- olmd <- score <- pathway <- labs <- V1 <- V2 <- labs <- cluster <- NULL
  
  
  set.seed(432143)
  
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
    }} else{
      df$score <- 1
      print("No valid gene set size columns provided. Please make sure input data has size_pathway or K, or else provide p-values to size points")
      scr_str <- "Gene Set"
    }
  
  if(is.null(cluster_result$cluster_membership$sign)){# if there is no sign column, result is from hypergeometric there is one cluster input
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
    
    plot <- local({
      . <- score <- pathway <- labs <- V1 <- V2 <- labs <- cluster <- NULL
      plot_df <- plot_df
      xlab <- xlab
      ylab <- ylab
      
      ggplot2::ggplot(plot_df) +
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
        ggplot2::scale_size_continuous(range = c(1, 10))
    })
    
    
  } else{
    # if there is a sign column, result is from GSEA and clustering is sign-separated
    figlist <- list()
    for(s in unique(cluster_result$cluster_membership$sign)){
      cl_df <- cluster_result$cluster_membership %>% 
        dplyr::filter(sign == s)
      df_temp <- df %>% 
        dplyr::filter(pathway %in% cl_df$pathway) %>% 
        dplyr::left_join(cl_df, by = c("pathway" = "pathway"))
      
      olmd <- cluster_result[[paste0(s, "_dist_mat")]]
      
      ## run dimension reduction on 2 axes ##
      if(dimred == "UMAP"){
        custom.settings = umap::umap.defaults
        custom.settings$input = "dist"
        custom.settings$random_state = 432143
        plot_df <- umap::umap(olmd, custom.settings)
        plot_df <- plot_df$layout %>% 
          as.data.frame() %>%
          tibble::rownames_to_column(var = "pathway") %>% 
          dplyr::left_join(df_temp, by = c("pathway" = "pathway"))
        
        xlab <- "UMAP1"
        ylab <- "UMAP2"
        
      } else if(dimred == "PCoA"){
        plot_df <- as.data.frame(stats::cmdscale(olmd)) %>% 
          tibble::rownames_to_column(var = "pathway") %>% 
          dplyr::left_join(df_temp, by = c("pathway" = "pathway"))
        
        xlab <- "PCoA1"
        ylab <- "PCoA2"
        
        
      } else if(dimred =="tSNE"){
        plot_df <- Rtsne::Rtsne(olmd, is_distance = TRUE, dims = 2, verbose = F, maxiter = 100)
        if(is.null(plot_df)){stop("Failed tSNE checks. Try a different dimred method")}
        plot_df <- as.data.frame(plot_df$Y) %>% 
          dplyr::mutate("pathway" = rownames(olmd), .before = "V1")%>% 
          dplyr::left_join(df_temp, by = c("pathway" = "pathway"))
        
        xlab <- "tSNE1"
        ylab <- "tSNE2"
        
      } else(stop("Valid options for dimred are 'UMAP', 'PCoA', or 'tSNE'."))
      
      ## make plot ##
      plot_df<- plot_df %>% 
        dplyr::mutate(score_scale = scales::rescale(score, to = c(1, 30)),
                      labs = stringr::str_replace_all(pathway, "_", " "),
                      labs = tolower(labs),
                      labs = stringr::str_wrap(labs, width = 25)) 
      
      figlist[[s]] <- local({
        . <- score <- pathway <- labs <- V1 <- V2 <- labs <- cluster <- NULL
        plot_df <- plot_df
        s <- s
        xlab <- xlab
        ylab <- ylab
        ggplot2::ggplot(plot_df) +
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
          ggplot2::scale_size_continuous(range = c(1, 10)) + 
          ggplot2::ggtitle(s)
      })
      
    }
    plot <- figlist %>% patchwork::wrap_plots(ncol = 1)
    
  }
  
  return(plot)
  
}