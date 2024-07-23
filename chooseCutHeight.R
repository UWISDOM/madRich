# df: result of enrichment function. Should contain a column of pathway names ("pathway"), p-values ("pvalue"), gene set category ("gs_cat"), and gene set subcategory where applicable ("gs_subcat"). If filtering by FDR, FDR values must be included as "FDR".
# category: a vector of Broad gene set categories included among the pathway names in df
# subcategory: a named vector of Broad gene set subcategories (names are corresponding cats provided in 'category')
# db: optional, custom gene set database formatted like MSigDB
# ID: "SYMBOL", "ENSEMBL", or "ENTREZ". What format of gene ID do you want to use for clustering.
# species: "human" or "mouse" for msigDB
# hclust_heights: vector of heights to test  for cuttree. Values must be between 0 and 1
# kmax: max number of clusters to assess for silhouette plot
# fdr_cutoff: optional, filter to apply to df before testing clustering solutions

chooseCutHeight <- function(
    df = NULL,
    category = NULL,
    subcategory = NULL,
    db = NULL,
    ID = "SYMBOL",
    species = "human",
    hclust_heights = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    kmax = NULL,
    fdr_cutoff = NULL
    
){
  db_list <- database <- db_format <- subcat <- NULL
  
  ### format df ###
  if(!is.null(fdr_cutoff)){
    df <- df %>% 
      dplyr::filter(FDR < fdr_cutoff)
  }
  
  ### get gene set matrix ###
  pathways <- unique(df$pathway)
  
  if(is.null(db)){
    db_list <- list()
    for(cat in category){
      database <- msigdbr::msigdbr(species, cat)
      if(cat == "C2"){
        database <- database %>% # subcategory in C2 is formatted bad. separate "CP" from "KEGG", etc.
          tidyr::separate(gs_subcat, sep = ":", into = c("gs_subcat_format", "x"), fill = "right")
      } else{
        database <- database %>%
          dplyr::mutate(gs_subcat_format = gs_subcat)
      }
      
      if(cat %in% names(subcategory)){ # filter to subcat when applicable
        database <- database %>% 
          dplyr::filter(gs_subcat_format == subcategory[[cat]])
      }
      
      
      database <- database %>% 
        dplyr::filter(gs_name %in% pathways) # remove pathways not in enrichment result
      db_list[[cat]] <- database
    }
    
    db_format <- dplyr::bind_rows(db_list)
    db_format$pathway <- db_format$gs_name
    
    if(ID == "SYMBOL") {
      db_format$gene <- db_format$gene_symbol
    } else if(ID == "ENSEMBL"){
      db_format$gene <- db_format$ensembl_gene
    } else if(ID == "ENTREZ"){
      db_format$gene <- db_format$entrez_gene
    } else{stop("ID error")}
    
  } else{
    db_format <- db
  }
  
  db_format <- db_format %>% 
    dplyr::select(c("pathway", "gene", "gs_description")) %>% 
    dplyr::distinct()
  
  
  ### make distance matrix ###
  dbl <- with(db_format, base::split(gene, pathway)) # convert pathway database to named list
  olm <- matrix(ncol = length(dbl), nrow = length(dbl)) # make nset x nset matrix
  colnames(olm) <- names(dbl)
  rownames(olm) <- names(dbl)
  
  for(i in names(dbl)){
    for(j in names(dbl)){
      if(i == j){
        v <- 1 # self gets 1
      }
      else{
        v <- length(base::intersect(dbl[[i]], dbl[[j]])) / min(c(length(dbl[[i]]), length(dbl[[j]]))) # calculate overlap coefficient
      }
      olm[i,j] <- v
    }
  }
  
  olmd <- 1-olm # convert similarity matrix to distance matrix
  
  ### make clusters ###
  cl <- stats::hclust(d = as.dist(olmd), method = "average")
  
  nclust <- c()
  
  for(i in hclust_heights){
    if(i > 1 | i < 0){stop("'hclust_heights' must be between 0 and 1.")}
    cl_cut <- stats::cutree(cl, h = i)
    cl_cut_df <- data.frame(cluster = as.factor(cl_cut))
    cl_cut_df$pathway <- names(cl_cut)
    nclust<- c(nclust, length(unique(cl_cut_df$cluster)))
  }
  
  clust_df <- data.frame("cut_height" = hclust_heights,
                         "k_at_height" = nclust)
  
  if(is.null(kmax)){
    kmax <- nrow(olm)-1
  }
  
  silhouette_score <- factoextra::fviz_nbclust(x = olm, diss = as.dist(olmd), FUN = factoextra::hcut, method = "silhouette", k.max=kmax)
  k_opt <- silhouette_score$data$clusters[which(silhouette_score$data$y == max(silhouette_score$data$y))]
  
  clust_df2 <- clust_df
  clust_df2[["y"]] <- (max(silhouette_score$data$y) - min(silhouette_score$data$y)) / 4
  silhouette_score <- silhouette_score +
    ggplot2::geom_vline(xintercept = clust_df2$k_at_height, linetype = "dashed") +
    ggplot2::geom_text(data = clust_df2, ggplot2::aes(x = k_at_height, y = y, label = paste0("height = ", cut_height, ", k = ", k_at_height)), angle = 90, vjust = 1.3) + 
    ggplot2::theme(text=ggplot2::element_text(size=10),
                  axis.text.x = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(color = "red"),
                  axis.ticks.x = ggplot2::element_blank()) +
    
    ggplot2::geom_vline(xintercept = k_opt, color = "red") +
    ggplot2::ggtitle("optimal: k = 22")
  
  clust_df3 <- clust_df %>% 
    dplyr::left_join(
      dplyr::mutate(silhouette_score$data, clusters = as.integer(clusters)), 
      by = c("k_at_height" = "clusters")) %>% 
    dplyr::rename("avg_silhouette_width" = y) %>%
    dplyr::arrange(cut_height)
  
  final <- list()
  final[["silhouette_score_plot"]] <- silhouette_score
  final[["summary_df"]] <- clust_df3
  final[["best_height"]] <- clust_df3$cut_height[which(clust_df3$avg_silhouette_width == max(clust_df3$avg_silhouette_width))]
  
  return(final)
}
