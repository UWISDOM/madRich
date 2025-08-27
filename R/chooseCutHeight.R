#' Plots to help guide choice of cut height for hierarchical clustering of gene sets
#'
#' @param df result of enrichment function. Should contain a column of pathway names ("pathway"), p-values ("pval"), gene set collections ("gs_collection"), and gene set subcollections where applicable ("gs_subcollection"). If filtering by FDR, k/K, or NES, these values must be included as "FDR", "k/K", or "NES", respectively.
#' @param enrich_method "hypergeometric" or "gsea." Serves as a sanity check for df filtering
#' @param collections a vector of Broad gene set collections included among the pathway names in df
#' @param subcollections a named vector of Broad gene set subcollections (names are corresponding cats provided in 'collections')
#' @param db optional, custom gene set database formatted like MSigDB
#' @param ID "SYMBOL", "ENSEMBL", or "ENTREZ". What format of gene ID do you want to use for clustering.
#' @param species "Homo sapiens" for human or "Mus musculus" for mouse
#' @param hclust_heights vector of heights to test  for cuttree. Values must be between 0 and 1
#' @param fdr_cutoff optional, FDR filter to apply to df before testing clustering solutions
#' @param abs_NES_cutoff optional, NES filter to apply to df before testing clustering solutions
#' @param kK_cutoff optional, k/K ratio filter to apply to df before testing clustering solutions
#' @param group_name optional, if your dataset contains multiple groups in a 'group' column, select pathways from provided group. Otherwise, all pathways passing filter will be clustered. 
#'
#' @return a list containing a plot(s) of silhouette scores by number of clusters, a summary dataframe, and a suggested "optimal" cut height
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
#' 
#' chooseCutHeight(df = df, enrich_method="hypergeometric",
#'                ID = "ENSEMBL",
#'                collections = c("C5"),
#'                subcollections = c("C5" = "GO:MF", "C5" = "GO:BP"),
#'                hclust_heights = c(0.3,0.5,0.7),
#'                group_name = "HRV1",
#'                fdr_cutoff = 0.4)

chooseCutHeight <- function(
    df = NULL,
    enrich_method = "hypergeometric",
    collections = NULL,
    subcollections = NULL,
    db = NULL,
    ID = "SYMBOL",
    species = "Homo sapiens",
    hclust_heights = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    fdr_cutoff = NULL, 
    abs_NES_cutoff = NULL, 
    kK_cutoff = NULL,
    group_name = NULL
    
){


. <- clusters <- cut_height <- y <- k_at_height <- gs_name <- gs_subcat_format <- gs_subcat <- group <- `k/K` <- NES <- FDR <- db_list <- database <- db_format <- subcat <- gs_subcollection <- gs_subcollection_format <- NULL
  

#Errors
if(!is.null(subcollections) & is.null(names(subcollections))){stop("subcollections must be a named vector.")}

  if(enrich_method == "gsea"){
    print("Note: when using for GSEA, this function will generate silhouette scores for all pass-filter gene sets separated by sign of NES. The 'best' solution will be the lower cut height of the two (e.g. err on the side of more clusters).")
  }

  ### format df ###
  if(!is.null(fdr_cutoff)){
    df <- df %>% 
      dplyr::filter(FDR < fdr_cutoff) # apply FDR cutoff
  }
  if(!is.null(abs_NES_cutoff)){ # apply NES cutoff if input is GSEA
    if(enrich_method != "gsea"){
      stop("NES filtering is only available with the 'gsea' method.")
    }
    df <- df %>% 
      dplyr::filter(abs(NES) > abs_NES_cutoff)

  }
  if(!is.null(kK_cutoff)){ # apply k/K cutoff if input is enrichment
    if(enrich_method != "hypergeometric"){
      stop("k/K filtering is only available with the 'hypergeometric' method.")
    }
    df <- df %>% 
      dplyr::filter(`k/K` > kK_cutoff)
  }
  if(!is.null(group_name)){ # filter by group of interest if provided
    df <- df%>% 
      dplyr::filter(group == group_name)
  }

  #check filtering
  if(length(unique(df$pathway)) == 0){stop("Looks like your enrichment result is empty after filtering. Adjust settings.")}
  
  
  ### get gene set matrix ###
  pw_list <- list()
  if(enrich_method == "hypergeometric") {
    pw_list[["pathways"]] <- unique(df$pathway)
  } else if(enrich_method == "gsea"){ # separate clustering by positive and negative part if GSEA
    pw_list[["positive"]] <- unique(df$pathway[which(df$NES > 0)])
    pw_list[["negative"]] <- unique(df$pathway[which(df$NES < 0)])
  } else{stop("Options for enrichment method are 'hypergeometric' and 'gsea'.")}
  
  
  if(is.null(db)){ # if database is not provided
    msigdbdf::msigdbdf("HS")
    db_list <- list()
    for(cat in collections){
      
      database <- msigdbr::msigdbr(species = species, collection = cat) # get msigdb reference
      if(cat == "C2"){
        database <- database %>% # subcollections in C2 is formatted bad. separate "CP" from "KEGG", etc.
          tidyr::separate_wider_delim(gs_subcollection, ":", names = c("gs_subcollection_format", "x"), too_few = "align_start", too_many = "merge")
      } else{
        database <- database %>%
          dplyr::mutate("gs_subcollection_format" = gs_subcollection)
      }
      
      if(cat %in% names(subcollections)){ # filter to subcat when applicable
      
      #Deal with multiple subcat per cat
      subcats <- subcollections[names(subcollections)==cat]
      if(length(unique(subcats))>1){
        # filter to subcat when applicable
        database <- database %>% 
          dplyr::filter(gs_subcollection_format == subcollections[[cat]])
      }
      }
      
      if(length(pw_list) > 1){# remove pathways not in enrichment result
        db_pos <- database %>% 
          dplyr::filter(gs_name %in% pw_list[["positive"]])%>% 
          dplyr::mutate(sign = "positive")
        db_neg <- database %>% 
          dplyr::filter(gs_name %in% pw_list[["negative"]])%>% 
          dplyr::mutate(sign = "negative")
        db_temp <- dplyr::bind_rows(db_pos, db_neg)
      } else if(length(pw_list) == 1){
        db_temp <- database %>% 
          dplyr::filter(gs_name %in% pw_list[["pathways"]]) %>% 
          dplyr::mutate(sign = "pathways")
      } else{stop()}
      db_list[[cat]] <- db_temp
    
    }
    
    db_format <- dplyr::bind_rows(db_list)
    
    if(ID == "SYMBOL") { # assign gene names based on indicated format
      db_format$gene <- db_format$gene_symbol
    } else if(ID == "ENSEMBL"){
      db_format$gene <- db_format$ensembl_gene
    } else if(ID == "ENTREZ"){
      db_format$gene <- db_format$entrez_gene
    } else{stop("ID error")}
    
  } else{ # if database is provided
    if(length(pw_list) > 1){# remove pathways not in enrichment result
      db_pos <- db %>% 
        dplyr::filter(gs_name %in% pw_list[["positive"]])%>% 
        dplyr::mutate(sign = "positive")
      db_neg <- db %>% 
        dplyr::filter(gs_name %in% pw_list[["negative"]])%>% 
        dplyr::mutate(sign = "negative")
      db_temp <- dplyr::bind_rows(db_pos, db_neg)
    } else if(length(pw_list) == 1){
      db_temp <- db %>% 
        dplyr::filter(gs_name %in% pw_list[["pathways"]]) %>% 
        dplyr::mutate(sign = "pathways")
    } else{stop()}
    db_format <- db_temp
  }
  db_format$pathway <- db_format$gs_name
  db_format <- db_format %>% 
    dplyr::select(c("sign", "pathway", "gene", "gs_description")) %>% 
    dplyr::distinct()
  
  final <- list()
  final[["database_format"]] <- db_format

  ### make distance matrix ###
  
  if(length(unique(db_format$sign)) > 1){
    plotlist <- list()
    dflist <- list()
    for(s in unique(db_format$sign)){ # separates if positive and negative signs
      db_temp <- db_format %>% 
        dplyr::filter(sign == s)
      dbl <- with(db_temp, base::split(gene, pathway)) # convert pathway database to named list
      olm <- matrix(ncol = length(dbl), nrow = length(dbl)) # make nset x nset matrix
      colnames(olm) <- names(dbl)
      rownames(olm) <- names(dbl)
      for(i in names(dbl)){ # iterate through matrix
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
      if(nrow(olmd)<=1){stop("Too few gene sets or no overlapping genes found across gene sets. Consider relaxing cutoffs to obtain more gene sets for clustering.")}
      cl <- stats::hclust(d = stats::as.dist(olmd), method = "average")
      
      nclust <- c()
      
      for(i in hclust_heights){ # cut tree for cut heights in input
        if(i > 1 | i < 0){stop("'hclust_heights' must be between 0 and 1.")}
        cl_cut <- stats::cutree(cl, h = i)
        cl_cut_df <- data.frame(cluster = as.factor(cl_cut))
        cl_cut_df$pathway <- names(cl_cut)
        nclust<- c(nclust, length(unique(cl_cut_df$cluster)))
      }
      

      clust_df <- data.frame("cut_height" = hclust_heights, # output df of n clusters vs cut height
                             "k_at_height" = nclust) 

      ## calculate cluster quality ##
      silhouette_score <- factoextra::fviz_nbclust(x = olm, diss = stats::as.dist(olmd), FUNcluster = factoextra::hcut, method = "silhouette", k.max= (nrow(olm) - 1))      
      k_opt <- silhouette_score$data$clusters[which(silhouette_score$data$y == max(silhouette_score$data$y))]
      clust_df2 <- clust_df
      clust_df2[["y"]] <- (max(silhouette_score$data$y) - min(silhouette_score$data$y)) / 4
      
      ## silhouette score plot ##
      plotlist[[s]] <- local({ 
        silhouette_score <- silhouette_score
        clust_df2 <- clust_df2
        k_opt <- k_opt
        silhouette_score +
          ggplot2::geom_vline(xintercept = clust_df2$k_at_height, linetype = "dashed") +
          ggplot2::geom_text(data = clust_df2, ggplot2::aes(x = k_at_height, y = y, label = paste0("height = ", cut_height, ", k = ", k_at_height)), angle = 90, vjust = 1.3) + 
          ggplot2::theme(text=ggplot2::element_text(size=10),
                         axis.text.x = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(color = "red"),
                         axis.ticks.x = ggplot2::element_blank()) +
          
          ggplot2::geom_vline(xintercept = k_opt, color = "red") +
          ggplot2::ggtitle(paste0(s, " NES optimal: k = ", k_opt))
      })
          
      clust_df3 <- clust_df %>% 
        dplyr::left_join(
          dplyr::mutate(silhouette_score$data, clusters = as.integer(clusters)), 
          by = c("k_at_height" = "clusters")) %>% 
        dplyr::rename("avg_silhouette_width" = y) %>%
        dplyr::arrange(cut_height)
      clust_df3$sign <- s
      dflist[[s]] <- clust_df3
      

    } 
    plot <- patchwork::wrap_plots(plotlist, ncol = 1)
    clust_df3 <- dplyr::bind_rows(dflist)
  
    
    } else{
      
    db_temp <- db_format 
    dbl <- with(db_temp, base::split(gene, pathway)) # convert pathway database to named list
    
    if(length(dbl) < 2){stop("Too few unique sets to cluster.")}
    
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
  if(nrow(olmd)<=1){stop("Too few gene sets or no overlapping genes found across gene sets. Consider relaxing cutoffs to obtain more gene sets for clustering.")}
  cl <- stats::hclust(d = stats::as.dist(olmd), method = "average")
  
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
  
  
  silhouette_score <- factoextra::fviz_nbclust(x = olm, diss = stats::as.dist(olmd), FUN = factoextra::hcut, method = "silhouette", k.max=(nrow(olm)-1))
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
    ggplot2::ggtitle(paste0("optimal: k = ", k_opt))
  
  clust_df3 <- clust_df %>% 
    dplyr::left_join(
      dplyr::mutate(silhouette_score$data, clusters = as.integer(clusters)), 
      by = c("k_at_height" = "clusters")) %>% 
    dplyr::rename("avg_silhouette_width" = y) %>%
    dplyr::arrange(cut_height)
  clust_df3$sign <- NA
  
  plot <- silhouette_score
  
  
  }
  
  final[["silhouette_score_plot"]] <- plot
  final[["summary_df"]] <- clust_df3
  final[["best_height"]] <- min(clust_df3$cut_height[which(clust_df3$avg_silhouette_width == max(clust_df3$avg_silhouette_width))])
  
  return(final)
}

