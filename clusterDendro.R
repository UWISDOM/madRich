#cluster_result: the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted

clusterDendro <- function(
    cluster_result = NULL
){
  
  set.seed(432143)
  
  dm <- cluster_result$dist_mat
  rownames(dm) <- stringr::str_replace_all(rownames(dm), "_", " ")
  rownames(dm) <- tolower(rownames(dm))
  colnames(dm) <- stringr::str_replace_all(colnames(dm), "_", " ")
  colnames(dm) <- tolower(colnames(dm))
  
  cl <- stats::hclust(d = stats::as.dist(dm), method = "average")
  
  dhc <- stats::as.dendrogram(cl)

  jdf <- cluster_result$cluster_membership %>% 
    mutate(label = stringr::str_replace_all(pathway, "_", " "),
           label = tolower(label),
           cluster = factor(cluster, levels = c(1:length(unique(cluster_result$cluster_membership$cluster)))))
  
  ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
  ddata$labels <- ddata$labels %>%
    dplyr::left_join(jdf, by = c("label" = "label")) %>% 
    dplyr::mutate(label = stringr::str_replace_all(label, "_", " "),
                  label = tolower(label))

  label_data <- dplyr::bind_cols(filter(ggdendro::segment(ddata), x == xend & x%%1 == 0), "label" = ddata$labels$label, "cluster" = ddata$labels$cluster) 

  dgdata <- ggdendro::segment(ddata) %>% 
    left_join(select(label_data, c(x, cluster))) %>% 
    mutate(lwd = NA)
  
    for(i in unique(jdf$cluster)){
      max <- max(dgdata$x[which(dgdata$cluster == i)])
      min <- min(dgdata$x[which(dgdata$cluster == i)])
      dgdata <- dgdata %>% 
        mutate(cluster = ifelse(x < max & x > min & xend < max & xend > min, i, cluster))
    }

  dgdata <- dgdata %>% 
    mutate(
      cluster = factor(cluster, levels = c(1:length(unique(cluster_result$cluster_membership$cluster)))))
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data=dgdata, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),lineend =  "round") +
    geom_text(data=label_data, y=0.04, hjust = 0, ggplot2::aes(x=xend, label=label), size=2) +
    coord_flip() + 
    scale_y_reverse(limits = c(NA,-0.5)) +
    ggplot2::geom_point(data = ddata$labels, ggplot2::aes(color = cluster, x = x), y = 0.02, size =4, alpha = 0.5, show.legend = F) +
    ggplot2::geom_text(data = ddata$labels, ggplot2::aes(label = cluster, x = x, color = cluster), y = 0.02, show.legend = F, size = 2)+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(), 
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()) 

  return(p)
}