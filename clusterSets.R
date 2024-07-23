# df: result of enrichment function. Should contain a column of pathway names ("pathway"), p-values ("pvalue"), gene set category ("gs_cat"), and gene set subcategory where applicable ("gs_subcat"). If filtering by FDR, FDR values must be included as "FDR".
# category: a vector of Broad gene set categories included among the pathway names in df
# subcategory: a named vector of Broad gene set subcategories (names are corresponding cats provided in 'category')
# db: optional, custom gene set database formatted like MSigDB
# ID: "SYMBOL", "ENSEMBL", or "ENTREZ". What format of gene ID do you want to use for clustering.
# species: "human" or "mouse" for msigDB
# hclust_height: Height for cutting tree for hclust. Value must be between 0 and 1
# fdr_cutoff: optional, filter to apply to df before clustering

clusterSets <- function(
    df = NULL,
    category = NULL,
    subcategory = NULL,
    db = NULL,
    ID = "SYMBOL",
    species = "human",
    hclust_height = 0.7,
    fdr_cutoff = NULL

){
  db_list <- database <- db_format <- subcat <- NULL
  
  final <- list()
  final[["input_df"]] <- df
  
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
  final[["datbase_format"]] <- db_format
  return(final)
}
