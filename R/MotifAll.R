#' Analysis of motifs in protein sequences
#'
#' This package implements the Motif_All algorithm
#'
#' @param x a vector of sequences
#' @param support the minimum support of motifs
#' @param k_min the minimum size of motifs
#'
#' @return a vector of motifs
#'
#' @export
#'
#' @author Guillaume Voisinne
#'
#' @examples
#' 
#' #load data :
#' x <- as.character(df_motif$Sequence[!is.na(df_motif$Cluster) & df_motif$Residue == "Y"])
#' motifs <- MotifAll(x, support = 0.05, k_min = 4)
#' 
MotifAll <- function(x, support=0.05, k_min=4){
  # x <- as.character(Tsite$Sequence[!is.na(Tsite$Cluster) & Tsite$Residue == "Y"])
  # 
  n_letters <- nchar(x[1])
  motif_list <- find_frequent_motifs(x = x, support = support)

  motifs <- unlist( lapply(motif_list, FUN = function(x){ write_motif(x, n_letters = n_letters) } ) )
  length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  unique_motifs <- unique(motifs[length_motifs > k_min])
  
}

#' @export
find_frequent_motifs <- function(x, support = 0.05) {
  # implementation of the Motif_All algorithm to find all frequent motifs (with a support >= k) in a set of sequences x
  n <- nchar(x[1])
  n_seq <- length(x)
  N_support <- dim(df)[1]*support
  
  M <- matrix("", n_seq, n)
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  
  df_bis<-as.data.frame(M, stringsAsFactors = FALSE)
  df<-as.data.frame(M)
  
  F_list <- list()
  
  # list all 1-motif
  count <- 0
  for (j in 1:n){
    sum_df <- summary(df[[j]])
    idx <- which( sum_df >= N_support )
    if(length(idx)>0){
      for (k in 1:length(idx)){
        count <- count + 1
        F_list[[count]] <- data.frame(positions = j, letters = names(sum_df)[idx[k]], stringsAsFactors = FALSE)
      }
    }
  }
  
  F_list_new <- F_list
  
  while ( length(F_list_new) > 0){
    
    count<-0
    F_list_1<- list()
    
    for (i in 1:length(F_list_new)){
      
      idx_match <- match_motif(df, F_list_new[[i]])
      
      for (j in setdiff(1:n, F_list_new[[i]]$positions)){
        sum_df <- summary( df[idx_match, j])
        idx <- which( sum_df >= N_support )
        if(length(idx)>0){
          for (k in 1:length(idx)){
            count <- count + 1
            F_list_1[[count]] <- rbind(F_list_new[[i]], 
                                       data.frame(positions = j, letters = names(sum_df)[idx[k]], stringsAsFactors = FALSE))
          }
        }
        
      }
      
    }
    
    F_list_new <- F_list_1
    F_list <- c(F_list, F_list_1)
  }
  
  return(F_list)
}

#' @export
match_motif <- function(df, motif){
  # get the indexes that match a given motif  
  df_int <- NULL
  for(i in length(motif$positions)){
    df_int <- cbind(df_int, df[ , motif$positions[i]] ==  motif$letters[i])
  }
  idx<-which(rowMeans(df_int)==1)
  return(idx)
  
}

#' @export
write_motif <- function(motif, n_letters){
  # write a motif as a string
  
  s<-rep(".", n_letters)
  for ( i in 1:length(motif$positions)){
    s[motif$positions[i]] <- motif$letters[i]
  }
  
  return(paste(s, collapse=""))
  
}
