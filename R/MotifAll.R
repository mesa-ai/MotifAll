#' Analysis of motifs in protein sequences
#'
#' This package implements the Motif_All algorithm
#'
#' @param x a vector of sequences (with a fixed length). Defined as the background set
#' @param idx_select Indexes of foreground sequences in \code{x}. Note that the 
#' foreground set must be a subset of the background set.
#' @param N_support minimum support of motifs in the foreground set (number of counts). 
#' Has priority over \code{support}
#' @param support minimum support of motifs in the foreground set (frequence)
#' @param signif significance level of output motifs (p-value associated with the Z-score)
#' @param var_p_val choose between p_value of Z-score ("p_value_Z_score")or p-value from hypergeometric test ("p_value_hyper")
#' @param k_min minimum size of motifs
#' @param k_max maximum size of motifs (set to NULL to ignore)
#' @param central_letter search only for motifs that have the given character in central position
#' @param match_central_letter Compute enrichment score on the restricted set of sequences that have the same central letter as the motif
#'
#' @return A list including the following elements :
#
#' @return \code{score} : data frame with the motifs identified along with scores and counts
#' @return \code{idx_match} : list of vectors containing the matching indexes of a given motif in the set of sequences x
#' 
#' @import stats
#' @import utils
#' 
#' @export
#'
#' @author Guillaume Voisinne
#'
#' @examples
#' #load data :
#' dir<- system.file("extdata", package = "MotifAll")
#' path <- paste(dir, "/phospho_site_annotated_plus_motifs.txt", sep="")
#' df_motif <- read.csv(path, sep="\t")
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
                     var_p_val  = "p_value_Z_score",
                     k_min = 3, 
                     k_max = NULL, 
                     central_letter = NULL,
                     match_central_letter=TRUE){
  # x : set of sequences
  # idx_select : indexes of foreground sequences in x
  
  n_letters <- nchar(x[1])
  motif_list <- find_frequent_motifs(x = x[idx_select], 
                                     N_support = N_support,
                                     support = support, 
                                     k_max = k_max, 
                                     central_letter = central_letter)
  motif_list <- filter_motif_list(motif_list, k_min = k_min)
  motifs <- print_motif(motif_list)
  #length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  #unique_motifs <- motifs[length_motifs >= k_min]
  
  score <- get_motif_score(x, idx_select, motif_list, match_central_letter = match_central_letter)
  idx_filter <- which(score[[var_p_val]] <= signif)
  idx_match <- vector("list", length(idx_filter))
  
  if(length(idx_filter)>0){
    
    score <- score[idx_filter, ]
    motif_list <- motif_list[idx_filter]
    idx_match <- match_motif_list_to_sequences(x, motif_list)
    
  } else {
    warning("No motif found. Try changing parameters.")
    return(NULL)
  }
  
  output <- list(score = score, motif_list = motif_list, idx_match = idx_match)
  
  return(output)
  
}

#' Converts a set of sequences to a dataframe
#'
#' @param x a vector of sequences (with a fixed length)
#' 
#' @return a data frame with letter positions as columns and sequences as rows
#' @export
sequence_to_df <- function(x){
  
  n <- nchar(x[1])
  n_seq<- length(x)
  
  M <- matrix("", n_seq, n)
  
  for (i in 1:n_seq){
    for (j in 1:n){
      M[i, j] <- substr(x[i],j,j)
    }
  }
  df <- as.data.frame(M)
  
  return(df)
}

#' Compute the enrichment score for a motif in the foreground set as compared to the background set
#'
#' @param df a data frame corresponding to set of background sequences as obtained by calling \code{sequence_to_df}. 
#' Could also be a vector of background sequences (with a fixed length) but will be slower.
#' @param idx_select Indexes of foreground sequences in \code{x}. Note that the foreground set must be a subset of the background set.
#' @param motif a motif
#' @param match_central_letter keep only sequences that have the same central letter as the motif
#' 
#' @return A data frame summarizing motif enrichment analysis in foreground set
#
#' @export
get_unique_motif_score <- function(df, idx_select, motif, match_central_letter=TRUE){
  
    if(!is.data.frame(df) & is.character(df)){
      df_int <- sequence_to_df(df) 
    } else {
      df_int <- df
    }
    
    motif_string <- print_motif(motif)
    n_center <- round(0.5 * (motif$size + 1))
    motif_central_letter <-  substr(motif_string, n_center, n_center)
    idx_match_center <- which(df_int[ ,n_center] == motif_central_letter)
      
    if(match_central_letter){
      idx_select_int <- intersect(idx_select, idx_match_center)
      df_select_int <- df_int[idx_select_int, ]
      df_int <- df_int[idx_match_center, ]
    } else {
      idx_select_int <- idx_select
      df_select_int <- df_int[idx_select_int, ]
    }
    
    n_sample <- length(idx_select_int)
    n_bckg <- dim(df_int)[1]

    motif_length<- length(motif$positions)
    
    idx_match_sample <- match_motif_df(df_select_int, motif)
    idx_match_bckg <-  match_motif_df(df_int, motif) 
    
    n_hits_sample <- length( idx_match_sample )
    n_hits_bckg <- length( idx_match_bckg )
    freq_sample <- n_hits_sample / n_sample
    freq_bckg <- n_hits_bckg / n_bckg
    
    c00 <- n_hits_sample
    c01 <- n_sample - c00
    c10 <- n_hits_bckg
    c11 <- n_bckg - c10
    
    LOR <- log(c00*c11/(c10*c01))
    SE <- sqrt(1/c00 + 1/c01 + 1/c10 + 1/c11)
    Z_score <-  LOR/SE
    p_value_Z_score <- 1-pnorm(Z_score)
    p_value_hyper <- 1-phyper(n_hits_sample-1, 
                                 n_hits_bckg,  
                                 n_bckg-n_hits_bckg,  
                                 n_sample);
    
    fold_change <- freq_sample/freq_bckg;
    
    score <- data.frame(motif = print_motif(motif),
                length = motif_length,
                Z_score = Z_score,
                p_value_Z_score = p_value_Z_score,
                fold_change = fold_change,
                p_value_hyper = p_value_hyper,
                n_hits_sample = n_hits_sample,
                n_sample = n_sample,
                freq_sample = freq_sample,
                n_hits_bckg = n_hits_bckg,
                n_bckg  = n_bckg,
                freq_bckg = freq_bckg)

    #output =  list(score=score, idx_match_bckg = idx_match_bckg, idx_match_sample = idx_match_sample)                
    
    return(score)
  
}

#' Compute the enrichment score for a list of motifs in the foreground set as compared to the background set
#'
#' @param x a vector of sequences (with a fixed length). Defined as the background set
#' @param idx_select Indexes of foreground sequences in \code{x}. Note that the foreground set must be a subset of the background set.
#' @param motif_list list of motifs
#' @param showProgress show progress bar in console
#' @param match_central_letter keep only sequences that have the same central letter as the motif
#' 
#' @return A data frame summarizing motif enrichment analysis in foreground set
#
#' @export
get_motif_score <- function(x, idx_select, motif_list, showProgress = TRUE, match_central_letter = TRUE){
  
  if(length(motif_list)==0){
    return(NULL)
  }
  
  df <- sequence_to_df(x)
  
  if (showProgress){
    cat("Computing enrichment score for each motif\n")
    pb <- txtProgressBar(min = 0, max = length(motif_list), style = 3)
  } 
  
  score <- vector("list", length(motif_list))
  for (i in 1:length(motif_list)){
    if (showProgress) setTxtProgressBar(pb, i)
    score[[i]] <- get_unique_motif_score(df, idx_select, motif_list[[i]], match_central_letter = match_central_letter)
  }
  if (showProgress) close(pb)
  
  output <- do.call(rbind, score)
  
  return(output)
  
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
          F_list[[count]] <- motif(positions = j, letters = names(idx[k]), size = n)
        }
      }
    }
    close(pb)
  } else {
    if( typeof(central_letter) == "character"){
      F_list[[1]] <- motif(positions = round((n+1)/2), letters = central_letter, size = n)
    }else{
      stop("Wrong parameter type: 'central_letter' should a 'character'")
    }
  }
  
  
  F_list_new <- get_unique_motifs(F_list)
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
            F_list_1[[count]] <- motif( positions = c(F_list_new[[i]]$positions, j),
                                        letters = c(F_list_new[[i]]$letters, names(sum_df)[idx[k]]),
                                        size = n )
                                      
          }
        }
        
      }
      
    }
    
    close(pb)
    
    if(length(F_list_1)>0){
      F_list_new <- get_unique_motifs(F_list_1)
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
match_sequence_to_motifs <- function(x, motif_list){
  
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

#' Match a list of motifs to a set of sequences
#' 
#' @param motif_list a list of motifs
#' @param x a vector of sequences (with a fixed length).
#' @param showProgress show progress bar in console
#' 
#' @return list of with matching indexes
#' 
#' @export
#' 
match_motif_list_to_sequences <- function(x, motif_list, showProgress=TRUE){
  
  df <- sequence_to_df(x)
  idx_match <- vector("list", length(motif_list))
  if (showProgress){
    cat("Matching motifs to sequences\n")
    pb <- txtProgressBar(min = 0, max = length(motif_list), style = 3)
  } 
  
  for (i in 1:length(motif_list)) {
    if (showProgress) setTxtProgressBar(pb, i)
    idx_match[[i]] <- match_motif_df(df, motif_list[[i]])
  }
  if (showProgress) close(pb)
  
  return(idx_match)
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

#' Converts a motif to a character string
#' 
#' @param x a motif
#' 
#' @return a character string
#' 
#' @export
motif_to_string<- function(x){
  
  s<-rep(".", x$size)
  for ( i in 1:length(x$positions)){
    s[x$positions[i]] <- x$letters[i]
  }
  
  return(paste(s, collapse=""))
}

#' Print a list motif as a vector of characters
#' 
#' @param motif_list a list of motifs
#' 
#' @return a vector of characters
#' 
#' @export
print_motif <- function(motif_list){
  if(class(motif_list)=="list"){
    motifs <- unlist( lapply(motif_list, FUN = function(x){ motif_to_string(x) } ) )
  } else {
    motifs <- motif_to_string(motif_list)
  }
  return(motifs)
}

#' Get unique motifs from a list of motifs
#' 
#' @param motif_list a list of motifs
#' 
#' @return a list of motifs
#' 
#' @export
get_unique_motifs <- function(motif_list){

  motifs <- print_motif(motif_list)
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

#' Converts a chqracter to a motif
#' 
#' @param s a character
#' @param null_char the character for unspecified letters
#' 
#' @return a motif
#' 
#' @export
convert_to_motif <- function(s, null_char = "."){
  
  n_letters <- nchar(s)
  s0<-strsplit(s, split="")[[1]]
  idx <- which(s0!=null_char)
  
  m <- motif( positions = idx, letters = s0[idx], size = n_letters)
  
  return(m)
}

#' Constructor for objects of class "motif"
#' 
#' @param positions a vector of integers specifying positions of letters in the motif
#' @param letters a vector of single characters (letters)
#' @param size the total number of characters (including unspecified letters)
#' 
#' @return a motif
#' 
#' @export
motif <- function( positions = NULL, letters = NULL, size = 1){
  
  m <- list(positions=positions, letters=letters, size = size)
  class(m) <- "motif"
  return(m)
  
}