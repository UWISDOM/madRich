#' Create word cloud plots to help in annotation of gene set clusters
#'
#' @param cluster_result the result of clusterSets which contains the overlap coefficient distance matrix and the cluster membership of gene sets to be plotted
#' @param rmwords a vector of words to remove from wordclouds. a standard list of words is included (articles, common broad category names, etc), but if your results have unexpected emphasis of useless words, you can add more here. 
#'
#' @return a list of ggplot objects
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
#' clusterWordclouds(res, rmwords = c("defense", "immune"))

clusterWordclouds <- function(
    cluster_result = NULL,
    rmwords = NULL
){
  
  #NOTE: a lot of this was adapted from vissE package (https://www.bioconductor.org/packages/release/bioc/html/vissE.html)
  
  . <- pathway <- gs_description <- cluster <- wd <- angle <- NULL
  set.seed(432143)
  plotlist <- list()
  
  #descriptions from database
  badwords <- c( # words not to include in frequency calculation
    'biocarta',
    'car',
    'cells',
    'dn',
    'expression',
    'gcm',
    'gene',
    'genes',
    'gnf2',
    'gobp',
    'gomf',
    'gocc',
    'gse',
    'hallmark',
    'kegg',
    'module',
    'morf',
    'neighborhood',
    'pathway',
    'pid',
    'reactome',
    'up',
    'wp',
    'is',
    'the',
    'medicus',
    'genes',
    'encode',
    'encoding',
    'of',
    'ensemble',
    'and',
    'for', 
    'by', 
    'important',
    'defining',
    'which',
    'system', 
    'via', 
    'proteins', 
    'to',
    'in',
    'upon',
    'or',
    'at',
    'a',
    'as',
    'such',
    'into',
    'their',
    'any',
    'that',
    'have',
    'from', 
    'its', 
    'be',
    'with', 
    'an',
    'over',
    'with', 
    'whose',
    'are',
    'an'
  )
  
  if(!is.null(rmwords)){ # add user-sumitted words
    badwords <- c(badwords, rmwords)
  }
  
  # hypergeometric - one result
  if(is.null(cluster_result$cluster_membership$sign)){
    
    df <- cluster_result$input_df
    cl_df <- cluster_result$cluster_membership
    
    df2 <- df %>%
      dplyr::left_join(cluster_result$database_format %>% dplyr::select(c(pathway, gs_description)) %>% dplyr::distinct()) %>%  # grab gene set 
      dplyr::left_join(cl_df, by = c("pathway" = "pathway"))
    
    gsv <- df2 %>% # get cluster list ordered by number of gene sets in cluster. 
      dplyr::group_by(cluster) %>%
      dplyr::count() %>%
      dplyr::arrange(-n) %>%
      dplyr::pull(cluster)
    
    wmp_n <- list() # holder for gene set names words
    wmp_d <- list() # holder for gene set description words
    wmp_b <- list()
    
    for(i in gsv){ # loop through clusters
      if(is.na(i)){next}
      
      dftemp <- df2 %>% # a lot of reformatting strings to remove weird characters, capitalization, etc.
        dplyr::filter(cluster == i) %>%
        dplyr::mutate(
          pathway = tolower(pathway),
          pathway = stringr::str_replace_all(pathway, "_", " "), 
          pathway = stringr::str_replace_all(pathway, "-", " "),
          pathway = stringr::str_replace_all(pathway, "/", " "),
          pathway = stringr::str_replace_all(pathway, ":", " "),
          pathway = stringr::str_replace_all(pathway, ",", ""),
          pathway = stringr::str_replace_all(pathway, stringr::fixed("["), " "),
          pathway = stringr::str_replace_all(pathway, stringr::fixed("]"), " "),
          pathway = stringr::str_replace_all(pathway, "[^[:alnum:]]", " "),
          gs_description = stringr::str_replace(gs_description, " \\s*\\[[^\\)]+\\]", ""),
          gs_description = tolower(gs_description),             
          gs_description = stringr::str_replace_all(gs_description, "-", " "),
          gs_description = stringr::str_replace_all(gs_description, "/", " "),
          gs_description = stringr::str_replace_all(gs_description, ":", " "),
          gs_description = stringr::str_replace_all(gs_description, ",", ""),
          gs_description = stringr::str_replace_all(gs_description, stringr::fixed("["), " "),
          gs_description = stringr::str_replace_all(gs_description, stringr::fixed("]"), " "),
          gs_description = stringr::str_replace_all(gs_description, "[^[:alnum:]]", " "),
          gs_description = stringr::str_replace_all(gs_description, stringr::fixed("."), "")
        )
      
      for(j in c("names", "descriptions", "both")){ # get counts of word freq
        
        if(j == "names"){
          cvec <- dftemp %>% # convert names to vectors of single words
            dplyr::pull(pathway) %>%
            unlist() %>%
            paste(collapse = " ") %>%
            stringr::str_split_1(" ") 
        } 
        else if(j == "descriptions"){
          cvec <-  dftemp %>% # convert descriptions to vectors of single words
            dplyr::pull(gs_description) %>% 
            unlist() %>%
            paste(collapse = " ") %>%
            stringr::str_split_1(" ") 
        }
        else if(j == "both"){
          cvec <- c(dftemp %>% # convert descriptions to vectors of single words
                      dplyr::pull(gs_description) %>% 
                      unlist() %>%
                      paste(collapse = " ") %>%
                      stringr::str_split_1(" "),
                    dftemp %>% # convert names to vectors of single words
                      dplyr::pull(pathway) %>%
                      unlist() %>%
                      paste(collapse = " ") %>%
                      stringr::str_split_1(" "))
        }
        
        cvec <- cvec[-which(cvec %in% badwords)]
        cdf <- data.frame("wd" = unique(cvec))
        
        
        tv <- c()
        for(k in cdf$wd){
          n <- length(cvec[which(cvec == k)])
          tv <- c(tv, n)
        }
        cdf$n <- tv
        cdf$angle = sample(c(0, 90), nrow(cdf), replace = TRUE, prob = c(0.65, 0.35))
        rl <- 5
        title <- paste0("Cluster ", i, ", GS ", j)
        
        # make wordcloud plot
        p3 <- ggplot2::ggplot(cdf, ggplot2::aes(label = wd, size = n, color = n, angle = angle)) +
          ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, shape = 'circle', eccentricity = 0.65) +
          scico::scale_colour_scico(palette = 'acton', direction = -1) +
          ggplot2::scale_size_area(max_size = 4 / log10(1 + n))+
          ggplot2::theme(
            panel.border = ggplot2::element_rect(colour = 'black', fill = NA),
            panel.grid = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = ggplot2::rel(rl) * 1.1),
            axis.text = ggplot2::element_text(size = ggplot2::rel(rl)),
            plot.title = ggplot2::element_text(size = ggplot2::rel(rl) * 0.5),
            strip.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            strip.text = ggplot2::element_text(size = ggplot2::rel(rl)),
            legend.text = ggplot2::element_text(size = ggplot2::rel(rl)),
            legend.title = ggplot2::element_text(size = ggplot2::rel(rl), face = 'italic') 
            
          ) +  ggplot2::ggtitle(title)
        
        if(j == "names"){
          wmp_n[[i]] <- p3
        } else if(j == "descriptions"){
          wmp_d[[i]] <- p3
        } else if(j == "both"){
          wmp_b[[i]] <- p3
        } else{stop("Please indicate what words you want to plot. Options include any or all of c('names', 'descriptions', 'both')")}
      }
    }
    
    
    plotlist[["wordcloud_names"]] <- wmp_n
    plotlist[["wordcloud_descriptions"]] <- wmp_d
    plotlist[["wordcloud_both"]] <- wmp_b
  }
  
  # gsea - list of 2 plots
  else{
    for(s in unique(cluster_result$cluster_membership$sign)){
      
      cl_df <- cluster_result$cluster_membership %>% 
        dplyr::filter(sign == s)
      df <- cluster_result$input_df %>% 
        dplyr::filter(pathway %in% cl_df$pathway)
      
      df2 <- df %>%
        dplyr::left_join(cluster_result$database_format %>% dplyr::select(c(pathway, gs_description)) %>% dplyr::distinct()) %>%  # grab gene set 
        dplyr::left_join(cl_df, by = c("pathway" = "pathway"))
      
      gsv <- df2 %>% # get cluster list ordered by number of gene sets in cluster. 
        dplyr::group_by(cluster) %>%
        dplyr::count() %>%
        dplyr::arrange(-n) %>%
        dplyr::pull(cluster)
      
      wmp_n <- list() # holder for gene set names words
      wmp_d <- list() # holder for gene set description words
      wmp_b <- list()
      
      for(i in gsv){ # loop through clusters
        if(is.na(i)){next}
        
        dftemp <- df2 %>% # a lot of reformatting strings to remove weird characters, capitalization, etc.
          dplyr::filter(cluster == i) %>%
          dplyr::mutate(
            pathway = tolower(pathway),
            pathway = stringr::str_replace_all(pathway, "_", " "), 
            pathway = stringr::str_replace_all(pathway, "-", " "),
            pathway = stringr::str_replace_all(pathway, "/", " "),
            pathway = stringr::str_replace_all(pathway, ":", " "),
            pathway = stringr::str_replace_all(pathway, ",", ""),
            pathway = stringr::str_replace_all(pathway, stringr::fixed("["), " "),
            pathway = stringr::str_replace_all(pathway, stringr::fixed("]"), " "),
            pathway = stringr::str_replace_all(pathway, "[^[:alnum:]]", " "),
            gs_description = stringr::str_replace(gs_description, " \\s*\\[[^\\)]+\\]", ""),
            gs_description = tolower(gs_description),             
            gs_description = stringr::str_replace_all(gs_description, "-", " "),
            gs_description = stringr::str_replace_all(gs_description, "/", " "),
            gs_description = stringr::str_replace_all(gs_description, ":", " "),
            gs_description = stringr::str_replace_all(gs_description, ",", ""),
            gs_description = stringr::str_replace_all(gs_description, stringr::fixed("["), " "),
            gs_description = stringr::str_replace_all(gs_description, stringr::fixed("]"), " "),
            gs_description = stringr::str_replace_all(gs_description, "[^[:alnum:]]", " "),
            gs_description = stringr::str_replace_all(gs_description, stringr::fixed("."), "")
          )
        
        for(j in c("names", "descriptions", "both")){ # get counts of word freq
          
          if(j == "names"){
            cvec <- dftemp %>% # convert names to vectors of single words
              dplyr::pull(pathway) %>%
              unlist() %>%
              paste(collapse = " ") %>%
              stringr::str_split_1(" ") 
          } 
          else if(j == "descriptions"){
            cvec <-  dftemp %>% # convert descriptions to vectors of single words
              dplyr::pull(gs_description) %>% 
              unlist() %>%
              paste(collapse = " ") %>%
              stringr::str_split_1(" ") 
          }
          else if(j == "both"){
            cvec <- c(dftemp %>% # convert descriptions to vectors of single words
                        dplyr::pull(gs_description) %>% 
                        unlist() %>%
                        paste(collapse = " ") %>%
                        stringr::str_split_1(" "),
                      dftemp %>% # convert names to vectors of single words
                        dplyr::pull(pathway) %>%
                        unlist() %>%
                        paste(collapse = " ") %>%
                        stringr::str_split_1(" "))
          }
          
          cvec <- cvec[-which(cvec %in% badwords)]
          cdf <- data.frame("wd" = unique(cvec))
          
          
          tv <- c()
          for(k in cdf$wd){
            n <- length(cvec[which(cvec == k)])
            tv <- c(tv, n)
          }
          cdf$n <- tv
          cdf$angle = sample(c(0, 90), nrow(cdf), replace = TRUE, prob = c(0.65, 0.35))
          rl <- 5
          title <- paste0("Cluster ", i, ", GS ", j)
          
          # make wordcloud plot
          p3 <- ggplot2::ggplot(cdf, ggplot2::aes(label = wd, size = n, color = n, angle = angle)) +
            ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, shape = 'circle', eccentricity = 0.65) +
            scico::scale_colour_scico(palette = 'acton', direction = -1) +
            ggplot2::scale_size_area(max_size = 4 / log10(1 + n))+
            ggplot2::theme(
              panel.border = ggplot2::element_rect(colour = 'black', fill = NA),
              panel.grid = ggplot2::element_blank(),
              axis.title = ggplot2::element_text(size = ggplot2::rel(rl) * 1.1),
              axis.text = ggplot2::element_text(size = ggplot2::rel(rl)),
              plot.title = ggplot2::element_text(size = ggplot2::rel(rl) * 0.5),
              strip.background = ggplot2::element_rect(fill = NA, colour = 'black'),
              strip.text = ggplot2::element_text(size = ggplot2::rel(rl)),
              legend.text = ggplot2::element_text(size = ggplot2::rel(rl)),
              legend.title = ggplot2::element_text(size = ggplot2::rel(rl), face = 'italic') 
              
            ) +  ggplot2::ggtitle(title)
          
          if(j == "names"){
            wmp_n[[i]] <- p3
          } else if(j == "descriptions"){
            wmp_d[[i]] <- p3
          } else if(j == "both"){
            wmp_b[[i]] <- p3
          } else{stop("Please indicate what words you want to plot. Options include any or all of c('names', 'descriptions', 'both')")}
        }
      }
      
      plotlist[[s]] <- list()
      plotlist[[s]][["wordcloud_names"]] <- wmp_n
      plotlist[[s]][["wordcloud_descriptions"]] <- wmp_d
      plotlist[[s]][["wordcloud_both"]] <- wmp_b
      
    }
    
  }
  return(plotlist)
}