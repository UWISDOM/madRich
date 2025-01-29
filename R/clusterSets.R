# df: result of enrichment function. Should contain a column of pathway names ("pathway"), p-values ("pvalue"), gene set categories ("gs_cat"), and gene set subcategories where applicable ("gs_subcat"). If filtering by FDR, k/K, or NES, these values must be included as "FDR", "k/K", or "NES", respectively.
# enrich_method: "hypergeometric" or "gsea." "gsea" will separate result by positive/negative NES within your filtered df.
# categories: a vector of Broad gene set categories included among the pathway names in df
# subcategories: a named vector of Broad gene set subcategories (names are corresponding cats provided in 'categories')
# db: optional, custom gene set database formatted like MSigDB
# ID: "SYMBOL", "ENSEMBL", or "ENTREZ". What format of gene ID do you want to use for clustering.
# species: "human" or "mouse" for msigDB
# hclust_height: Height for cutting tree for hclust. Value must be between 0 and 1
# fdr_cutoff: optional, high-end filter to apply to df before clustering
# abs_NES_cutoff: optional, low-end. absolute value filter to apply to df before clustering of GSEA results
# kK_cutoff: optional, low-end filter to apply to df before clustering of hypergeometric results
# group_name: optional, if your dataset contains multiple groups in a 'group' column, select pathways from provided group. Otherwise, all pathways passing filter will be clustered. 

#' Perform Hierarchical Clustering of Gene Sets Based on the Overlap Coefficient
#'
#' @param df result of enrichment function. Should contain a column of pathway names ("pathway"), p-values ("pvalue"), gene set categories ("gs_cat"), and gene set subcategories where applicable ("gs_subcat"). If filtering by FDR, k/K, or NES, these values must be included as "FDR", "k/K", or "NES", respectively.
#' @param enrich_method "hypergeometric" or "gsea." "gsea" will separate result by positive/negative NES within your filtered df.
#' @param categories a vector of Broad gene set categories included among the pathway names in df
#' @param subcategories a named vector of Broad gene set subcategories (names are corresponding cats provided in 'categories')
#' @param db optional, custom gene set database formatted like MSigDB
#' @param ID "SYMBOL", "ENSEMBL", or "ENTREZ". What format of gene ID do you want to use for clustering.
#' @param species "human" or "mouse" for msigDB
#' @param hclust_height Height for cutting tree for hclust. Value must be between 0 and 1
#' @param fdr_cutoff optional, high-end filter to apply to df before clustering
#' @param abs_NES_cutoff optional, low-end. absolute value filter to apply to df before clustering of GSEA results
#' @param kK_cutoff optional, low-end filter to apply to df before clustering of hypergeometric results
#' @param group_name optional, if your dataset contains multiple groups in a 'group' column, select pathways from provided group. Otherwise, all pathways passing filter will be clustered. 
#'
#' @return a list containing the input data frame, formatted reference database, cluster membership, and distance matrix/matrices. To be used in downstream visualizations.
#' @export
#'
#' @examples
#' # Run GSEA using SEARchways
#' gene_list <- SEARchways::example.gene.list
#' df1 <- SEARchways::BIGsea(gene_list=gene_list, 
#'              category="C2", subcategory="CP", ID="ENSEMBL")
#' df2 <- SEARchways::BIGsea(gene_list=gene_list, 
#'              category="C5", subcategory="GO:BP", ID="ENSEMBL")
#' df <- dplyr::bind_rows(df1, df2)
#' 
#' res <- clusterSets(df = df, enrich_method="hypergeometric",
#'                 ID = "ENSEMBL",
#'                 categories = c("C2","C5"),
#'                 #subcategories = c("C2" = "CP", "C5" = "GO:BP"),
#'                 hclust_height = 0.3,
#'                 group_name = "HRV1")
#' 
#' # Run enrichment using SEARchways
#' gene_list2 <- list(HRV1 = names(SEARchways::example.gene.list[[1]]),
#'                   HRV2 = names(SEARchways::example.gene.list[[2]]))
#' df1 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'              category="C5", subcategory="GO:MF", ID="ENSEMBL")
#' df2 <- SEARchways::BIGprofiler(gene_list=gene_list2, 
#'              category="C5", subcategory="GO:BP", ID="ENSEMBL")
#' df <- dplyr::bind_rows(df1, df2)
#' 
#' res <- clusterSets(df = df, enrich_method="hypergeometric",
#'                 ID = "ENSEMBL",
#'                 categories = c("C5"),
#'                 subcategories = c("C5" = "GO:MF", "C5" = "GO:BP"),
#'                 hclust_height = 0.3,
#'                 group_name = "HRV1",
#'                 fdr_cutoff = 0.3)

clusterSets <- function(
    df = NULL,
    enrich_method = "hypergeometric",
    categories = NULL,
    subcategories = NULL,
    db = NULL,
    ID = "SYMBOL",
    species = "human",
    hclust_height = 0.7,
    fdr_cutoff = NULL,
    abs_NES_cutoff = NULL, 
    kK_cutoff = NULL,
    group_name = NULL
){
  . <- gs_name <- gs_subcat_format <- gs_subcat <- group <- `k/K` <- NES <- FDR <- db_list <- database <- db_format <- subcat <- NULL
  
  #Errors
  if(!is.null(subcategories) & is.null(names(subcategories))){stop("subcategories must be a named vector.")}
  
  final <- list()
  final[["input_df"]] <- df
  
  ### format df ###
  if(!is.null(fdr_cutoff)){
    df <- df %>% 
      dplyr::filter(FDR < fdr_cutoff)
  }
  if(!is.null(abs_NES_cutoff)){
    if(enrich_method != "gsea"){
      stop("NES filtering is only available with the 'gsea' method.")
      }
    df <- df %>% 
      dplyr::filter(abs(NES) > abs_NES_cutoff)
  }
  if(!is.null(kK_cutoff)){
    if(enrich_method != "hypergeometric"){
      stop("k/K filtering is only available with the 'hypergeometric' method.")
      }
    df <- df %>% 
      dplyr::filter(`k/K` > kK_cutoff)
  }
  if(!is.null(group_name)){
    df <- df%>% 
      dplyr::filter(group == group_name)
  }
  
  #check filtering
  if(length(unique(df$pathway)) == 0){stop("Looks like your enrichment result is empty after filtering. Adjust settings.")}
    
  ### get gene set matrix ###
  pw_list <- list()
  if(enrich_method == "hypergeometric") {
    pw_list[["pathways"]] <- unique(df$pathway)
  } else if(enrich_method == "gsea"){
    pw_list[["positive"]] <- unique(df$pathway[which(df$NES > 0)])
    pw_list[["negative"]] <- unique(df$pathway[which(df$NES < 0)])
  } else{stop("Options for enrichment method are 'hypergeometric' and 'gsea'.")}
 
  if(is.null(db)){
    db_list <- list()
    for(cat in categories){
      database <- msigdbr::msigdbr(species, cat)
      
      # Fix C2 subcat names
      if(cat == "C2"){
        database <- database %>% # subcategories in C2 is formatted bad. separate "CP" from "KEGG", etc.
          tidyr::separate(gs_subcat, sep = ":", into = c("gs_subcat_format", "x"), fill = "right")
      } else{
        database <- database %>%
          dplyr::mutate(gs_subcat_format = gs_subcat)
      }
      
      #Deal with multiple subcat per cat
      subcats <- subcategories[names(subcategories)==cat]
      if(length(subcats)>0){
        # filter to subcat when applicable
        database <- database %>% 
          dplyr::filter(gs_subcat_format %in% subcats)
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
   
    if(ID == "SYMBOL") {
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
  final[["datbase_format"]] <- db_format
  
  
  ### make distance matrix ###
  cluster_list <- list()
  if(length(unique(db_format$sign)) > 1){
    for(s in unique(db_format$sign)){
      db_temp <- db_format %>% 
        dplyr::filter(sign == s)
      dbl <- with(db_temp, base::split(gene, pathway)) # convert pathway database to named list
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
      
      if(hclust_height > 1 | hclust_height < 0){stop("'hclust_height' must be between 0 and 1.")}
      else{
        cl_cut <- stats::cutree(cl, h = hclust_height)
        cl_cut_df <- data.frame(cluster = as.factor(cl_cut))
        cl_cut_df$pathway <- names(cl_cut)
        rownames(cl_cut_df) <- NULL
      }
      
      cl_cut_df$sign <- s
      cluster_list[[s]] <- cl_cut_df
      final[[paste0(s,"_dist_mat")]] <- olmd

    }
    cl_df_join <- dplyr::bind_rows(cluster_list)
    final[["cluster_membership"]] <- cl_df_join
    
  } else{
    db_temp <- db_format
    dbl <- with(db_temp, base::split(gene, pathway)) # convert pathway database to named list
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
    
    if(hclust_height > 1 | hclust_height < 0){stop("'hclust_height' must be between 0 and 1.")}
    else{
      cl_cut <- stats::cutree(cl, h = hclust_height)
      cl_cut_df <- data.frame(cluster = as.factor(cl_cut))
      cl_cut_df$pathway <- names(cl_cut)
      rownames(cl_cut_df) <- NULL
    }
    final[["dist_mat"]] <- olmd
    final[["cluster_membership"]] <- cl_cut_df
  }

  
  


  return(final)
}
