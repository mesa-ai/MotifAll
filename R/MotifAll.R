#' Analysis of motifs in protein sequences
#'
#' This package implements the Motif_All algorithm
#'
#' @param x a vector of sequences (with a fixed length). Defined as the background set
#' @param idx_select Indexes of foreground sequences in \code{x}. Note that the foreground set must be a subset of the background set.
#' @param N_support minimum support of motifs in the foreground set (number of counts). Has priority over \code{support}
#' @param support minimum support of motifs in the foreground set (frequence)
#' @param signif significance level of output motifs (p-value associated with the Z-score)
#' @param k_min minimum size of motifs
#' @param k_max maximum size of motifs (set to NULL to ignore)
#' @param central_letter search only for motifs that have the given character in central position
#'
#' @return A list including the following elements :
#
#' @return \code{score} : data frame with the motifs identified along with scores and counts
#' @return \code{idx_match_bckg} : list of vectors containing the matching indexes of a given motif in the set of sequences x
#' @return \code{idx_match_sample} : list of vectors containing the matching indexes of a given motif in the set of sequences x also part of the foreground set
#' 
#' @import stats
#' @import utils
#' 
#' @export
#'
#' @author Guillaume Voisinne
#'
#' @examples
#' 
#' #load data :
#' dir<- system.file("extdata", package = "MotifAll")
#' path <- paste(dir,"/phospho_site_annotated_plus_motifs.txt",sep="")
#' df_motif <- read.csv(path, sep="\t", quote = "\"")
#' #data(df_motif, package = "MotifAll")
#' psite <- "S"
#' 
#' idx_residue <- which(substr(as.character(df_motif$Sequence),8,8) == psite)
#' df <- df_motif[idx_residue, ]
#' 
#' x <- as.character(df$Sequence)
#' idx_select <- which(!is.na(df$Cluster))
#' 
#' res <- MotifAll(x, idx_select, support = 0.05, signif = 1e-4, k_min = 1, central_letter = psite)
#' 
MotifAll <- function(x, 
                     idx_select,
                     N_support = NULL,
                     support = 0.05,
                     signif = 0.05,
                     k_min = 3, 
                     k_max = NULL, 
                     central_letter = NULL){
  # x : set of sequences
  # idx_select : indexes of foreground sequences in x
  
  n_letters <- nchar(x[1])
  motif_list <- find_frequent_motifs(x = x[idx_select], 
                                     N_support = N_support,
                                     support = support, 
                                     k_max = k_max, 
                                     central_letter = central_letter)
  motif_list <- filter_motif_list(motif_list, k_min = k_min)
  motifs <- print_motifs(motif_list, n_letters)
  #length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  #unique_motifs <- motifs[length_motifs >= k_min]
  
  res <- get_motif_score(x, idx_select, motif_list)
  idx_filter <- which(res$score$p_value_Z_score <= signif)
  if(length(idx_filter)>0){
    res$score <- res$score[idx_filter, ]
    res$idx_match_bckg <- res$idx_match_bckg[idx_filter]
    res$idx_match_sample <- res$idx_match_sample[idx_filter]
    res$motif_list <- motif_list[idx_filter]
  } else {
    warning("No motif found. Try changing parameters.")
    return(NULL)
  }
  
  return(res)
  
}


#' Compute the enrichment score for a list of motifs in the foreground set as compared to the background set
#'
#' @param x a vector of sequences (with a fixed length). Defined as the background set
#' @param idx_select Indexes of foreground sequences in \code{x}. Note that the foreground set must be a subset of the background set.
#' @param motif_list list of motifs
#'
#' @return A data frame summarizing motif enrichment analysis in foreground set
#
#' @export
get_motif_score <- function(x, idx_select, motif_list){
  
  if (length(motif_list)==0){
    warning("No motif found. Try changing parameters.")
    return(NULL)
  }
  
  n <- nchar(x[1])
  n_seq<- length(x)
  M <- matrix("", n_seq, n)
  
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  
  df <- as.data.frame(M)
  df_select <- df[idx_select, ]
  
  n_sample <- length(idx_select)
  n_bckg <- n_seq
  
  
  n_motifs <- length(motif_list)
  motif <- rep("", n_motifs)
  motif_length <- rep(NA, n_motifs)
  n_hits_sample <- rep(NA, n_motifs)
  n_hits_bckg <- rep(NA, n_motifs)
  freq_sample <- rep(NA, n_motifs)
  freq_bckg <- rep(NA, n_motifs)
  p_value_hyper <- rep(NA, n_motifs)
  fold_change <- rep(NA, n_motifs)
  Z_score<-  rep(NA, n_motifs)
  p_value_Z_score<- rep(NA, n_motifs)
  
  idx_match_sample <- list()
  idx_match_bckg <- list()
  
  cat("Computing motif scores\n")
  pb <- txtProgressBar(min = 0, max = length(motif_list), style = 3)
  
  for(i in 1:length(motif_list)){
    
    setTxtProgressBar(pb, i)
    motif[i] <- print_motifs(motif_list[i], n)
    motif_length[i] <- length(motif_list[[i]]$positions)
    idx_match_sample[[i]] <- idx_select[match_motif_df(df_select, motif_list[[i]])]
    idx_match_bckg[[i]] <- match_motif_df(df, motif_list[[i]])
    
     
    n_hits_sample[i] <- length( idx_match_sample[[i]] )
    n_hits_bckg[i] <- length( idx_match_bckg[[i]] )
    freq_sample[i] <- n_hits_sample[i] / n_sample
    freq_bckg[i] <- n_hits_bckg[i] / n_bckg
    
    c00 <- n_hits_sample[i]
    c01 <- n_sample - c00
    c10 <- n_hits_bckg[i]
    c11 <- n_bckg - c10
    
    LOR <- log(c00*c11/(c10*c01))
    SE <- sqrt(1/c00 + 1/c01 + 1/c10 + 1/c11)
    Z_score[i] <-  LOR/SE
    p_value_Z_score[i] <- 1-pnorm(Z_score[i])
    p_value_hyper[i] <- 1-phyper(n_hits_sample[i]-1, 
                          n_hits_bckg[i],  
                          n_bckg-n_hits_bckg[i],  
                          n_sample);
    
    fold_change[i] <- freq_sample[i]/freq_bckg[i];
    
  }
  close(pb)
  
  df <- data.frame(motif = motif,
                   length = motif_length,
                   Z_score = Z_score,
                   p_value_Z_score = p_value_Z_score,
                   fold_change = fold_change,
                   p_value_hyper = p_value_hyper,
                   n_hits_sample = n_hits_sample,
                   n_sample = rep(n_sample, n_motifs),
                   freq_sample = freq_sample,
                   n_hits_bckg = n_hits_bckg,
                   n_bckg  = rep(n_bckg, n_motifs),
                   freq_bckg = freq_bckg
                   )
  
  res <- list(score=df, idx_match_bckg = idx_match_bckg, idx_match_sample = idx_match_sample)
  
  return(res)
  
}

#' Filter motifs based on their length
#'
#' @param motif_list list of motifs
#' @param k_min minimum motif length
#'
#' @return A list of motifs
#
#' @export
filter_motif_list <- function(motif_list, k_min){

  length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  idx_filter = which(length_motifs >= k_min)
  motif_list_filter <- motif_list[idx_filter]

  return(motif_list_filter)
  
}

#' Find frequent motifs in a set of sequences using the Motif_All algorithm
#'
#' @param x a vector of sequences (with a fixed length). Defined as the background set
#' @param N_support minimum support of motifs in the foreground set (number of counts). Has priority over \code{support}
#' @param support minimum support of motifs in the foreground set (frequence)
#' @param k_max maximum size of motifs (set to NULL to ignore)
#' @param central_letter search only for motifs that have the given character in central position
#'
#' @return a list of motifs
#' 
#' @export
find_frequent_motifs <- function(x, N_support = NULL, support = 0.05, k_max = NULL, central_letter = NULL) {
  # implementation of the Motif_All algorithm to find all frequent motifs (with a support >= k) in a set of sequences x
  
  n <- nchar(x[1])
  
  if(!is.null(k_max)){
    k_lim <- k_max
  } else {
    k_lim <- n
  }
    
  n_seq <- length(x)
  
  if (is.null(N_support)){
    N_support_eff <- n_seq*support
  } else {
    N_support_eff <- N_support
  }
  
  
  M <- matrix("", n_seq, n)
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  
  df_bis<-as.data.frame(M, stringsAsFactors = FALSE)
  df<-as.data.frame(M)
  
  F_list <- list()
  
  if(is.null(central_letter)){
    
    # list all 1-motif
    count <- 0
    cat("Listing 1-motif\n")
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    
    for (j in 1:n){
      setTxtProgressBar(pb, j)
      sum_df <- summary(df[[j]])
      idx <- which( sum_df >= N_support_eff )
      if(length(idx)>0){
        for (k in 1:length(idx)){
          count <- count + 1
          F_list[[count]] <- data.frame(positions = j, letters = names(idx[k]), stringsAsFactors = FALSE)
        }
      }
    }
    close(pb)
  } else {
    if( typeof(central_letter) == "character"){
      F_list[[1]] <- data.frame(positions = round((n+1)/2), letters = central_letter, stringsAsFactors = FALSE)
    }else{
      stop("Wrong parameter type: 'central_letter' should a 'character'")
    }
  }
  
  
  F_list_new <- get_unique_motifs(F_list, n_letters = n)
  cat(paste("N=",length(F_list_new),"\n",sep=""))
  #F_list_new <- F_list
  
  k_motif <- 1
  while ( length(F_list_new) > 0 & k_motif < k_lim){
    k_motif <- k_motif + 1
    cat(paste("Listing ",k_motif,"-motif from ", length(F_list_new), " parent motifs\n", sep=""))
    pb <- txtProgressBar(min = 0, max = length(F_list_new), style = 3)
    
    count<-0
    F_list_1 <- list()
    
    for (i in 1:length(F_list_new)){
      setTxtProgressBar(pb, i)
      idx_match <- match_motif_df(df, F_list_new[[i]])
      
      for (j in setdiff(1:n, F_list_new[[i]]$positions)){
        sum_df <- summary( df[idx_match, j])
        idx <- which( sum_df >= N_support_eff )
        if(length(idx)>0){
          for (k in 1:length(idx)){
            count <- count + 1
            F_list_1[[count]] <- rbind(F_list_new[[i]], 
                                       data.frame(positions = j, letters = names(sum_df)[idx[k]], stringsAsFactors = FALSE))
          }
        }
        
      }
      
    }
    
    #print( print_motifs(F_list_1, n) )
    close(pb)
    
    if(length(F_list_1)>0){
      F_list_new <- get_unique_motifs(F_list_1, n_letters = n)
    } else {
      F_list_new <- list()
    }
    
    cat(paste("N=",length(F_list_new),"\n",sep=""))
    
    F_list <- c(F_list, F_list_new)
  }
  
  return(F_list)
}

#' Match a sequence to a list of motifs
#' 
#' @param motif_list list of motifs
#' @param x a sequence
#' 
#' @return vector of matching indexes
#' 
#' @export
match_sequence_to_motifs <- function(motif_list, x){
  
  n <- nchar(x)
  n_seq <- 1
  
  M <- matrix("", n_seq, n)
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  
  df<-as.data.frame(M)
  
  idx<-NULL
  for (i in 1:length(motif_list)){
    if (length(match_motif_df(df, motif_list[[i]])) > 0 ){
      idx <- c(idx, i)
    }
  }
  
  return(idx)
  
}

#' Match a motif to a set of sequences
#' 
#' @param motif a motif
#' @param x a vector of sequences (with a fixed length).
#' 
#' @return vector of matching indexes
#' 
#' @export
match_motif_to_sequences <- function(x, motif){
  
  n <- nchar(x[1])
  n_seq <- length(x)

  M <- matrix("", n_seq, n)
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  
  df<-as.data.frame(M)
  
  idx<-match_motif_df(df, motif)
  
  return(idx)
  
}

#' Match a motif to a set of sequences represented as a data frame
#' 
#' @param motif a motif
#' @param df a data frame corresponding to a set of sequences (with a fixed length).
#' 
#' @return vector of matching indexes
#' 
#' @export
match_motif_df <- function(df, motif){
  # get the indexes that match a given motif  
  df_int <- NULL
  for(i in 1:length(motif$positions)){
    df_int <- cbind(df_int, df[ , motif$positions[i]] ==  motif$letters[i])
  }
  idx<-which(rowMeans(df_int)==1)
  return(idx)
  
}

#' Print a motif as a character string
#' 
#' @param motif a motif
#' @param n_letters length of the output character string
#' 
#' @return a character string
#' 
#' @export
print_single_motif <- function(motif, n_letters){
  # write a motif as a string
  
  s<-rep(".", n_letters)
  for ( i in 1:length(motif$positions)){
    s[motif$positions[i]] <- motif$letters[i]
  }
  
  return(paste(s, collapse=""))
}

#' Print a list motif as a vector of characters
#' 
#' @param motif_list a list of motifs
#' @param n_letters length of the output character string
#' 
#' @return a vector of characters
#' 
#' @export
print_motifs <- function(motif_list, n_letters){
  motifs <- unlist( lapply(motif_list, FUN = function(x){ print_single_motif(x, n_letters = n_letters) } ) )
  return(motifs)
}

#' Get unique motifs from a list of motifs
#' 
#' @param motif_list a list of motifs
#' @param n_letters maximum number of letters in a motif
#' 
#' @return a list of motifs
#' 
#' @export
get_unique_motifs <- function(motif_list, n_letters){

  motifs <- print_motifs(motif_list, n_letters)
  unique_motifs <- unique(motifs)
  
  
  motif_list_unique <- list()
  #n_rep <- rep(NA, length(unique_motifs))
  #idx_rep <- list()
  for (i in 1:length(unique_motifs)){
    idx_u <- which(motifs == unique_motifs[i])
    #n_rep[i]=length(idx_u)
    #idx_rep[[i]] <- idx_u
    motif_list_unique[[i]] <- motif_list[[ idx_u[1] ]]
  }
  
  return(motif_list_unique)
}
