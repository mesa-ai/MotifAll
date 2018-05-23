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
#' data(df_motif, package = "MotifAll")
#' site_residue="S"
#' df <- df_motif[df_motif$Residue == site_residue, ]
#' 
#' x <- as.character(df$Sequence)
#' idx_select <- which(!is.na(df$Cluster))
#' 
#' res <- MotifAll(x, idx_select, support = 0.05, k_min = 3, central_letter = site_residue)
#' 
MotifAll <- function(x, 
                     idx_select, 
                     support=0.05, 
                     k_min=3, 
                     k_max = NULL, 
                     central_letter = NULL){
  # x : set of sequences
  # idx_select : indexes of foreground sequences in x
  
  n_letters <- nchar(x[1])
  motif_list <- find_frequent_motifs(x = x[idx_select], 
                                     support = support, 
                                     k_max = k_max, 
                                     central_letter = central_letter)
  motif_list <- filter_motif_list(motif_list, k_min = k_min)
  motifs <- print_motifs(motif_list, n_letters)
  #length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  #unique_motifs <- motifs[length_motifs >= k_min]
  
  res <- get_motif_score(x, idx_select, motif_list)
  
  return(res)
  
}

#' @export
get_motif_score <- function(x, idx_select, motif_list){
  
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
  p_value <- rep(NA, n_motifs)
  fold_change <- rep(NA, n_motifs)
  
  idx_match_sample <- list()
  idx_match_bckg <- list()
  
  for(i in 1:length(motif_list)){
    
    motif[i] <- print_motifs(motif_list[i], n)
    motif_length[i] <- length(motif_list[[i]]$positions)
    idx_match_sample[[i]] <- idx_select[match_motif(df_select, motif_list[[i]])]
    idx_match_bckg[[i]] <- match_motif(df, motif_list[[i]])
    n_hits_sample[i] <- length( idx_match_sample[[i]] )
    n_hits_bckg[i] <- length( idx_match_bckg[[i]] )
    freq_sample[i] <- n_hits_sample[i] / n_sample
    freq_bckg[i] <- n_hits_bckg[i] / n_bckg
    
    p_value[i] = 1-phyper(n_hits_sample[i]-1, 
                          n_hits_bckg[i],  
                          n_bckg-n_hits_bckg[i],  
                          n_sample);
    
    fold_change[i] = freq_sample[i]/freq_bckg[i];
    
  }
  
  df <- data.frame(motif = motif,
                   length = motif_length,
                   p_value = p_value,
                   fold_change = fold_change,
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

#' @export
filter_motif_list <- function(motif_list, k_min){

  length_motifs <- unlist( lapply(motif_list, FUN = function(x){ length(x$positions) } ) ) 
  idx_filter = which(length_motifs >= k_min)
  motif_list_filter <- motif_list[idx_filter]

  return(motif_list_filter)
  
}
                                 
#' @export
find_frequent_motifs <- function(x, support = 0.05, k_max = NULL, central_letter = NULL) {
  # implementation of the Motif_All algorithm to find all frequent motifs (with a support >= k) in a set of sequences x
  
  
  n <- nchar(x[1])
  
  if(!is.null(k_max)){
    k_lim <- k_max
  } else {
    k_lim <- n
  }
    
  n_seq <- length(x)
  N_support <- n_seq*support
  
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
      idx <- which( sum_df >= N_support )
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

#' @export
match_motif <- function(df, motif){
  # get the indexes that match a given motif  
  df_int <- NULL
  for(i in 1:length(motif$positions)){
    df_int <- cbind(df_int, df[ , motif$positions[i]] ==  motif$letters[i])
  }
  idx<-which(rowMeans(df_int)==1)
  return(idx)
  
}

#' @export
print_single_motif <- function(motif, n_letters){
  # write a motif as a string
  
  s<-rep(".", n_letters)
  for ( i in 1:length(motif$positions)){
    s[motif$positions[i]] <- motif$letters[i]
  }
  
  return(paste(s, collapse=""))
}

#' @export
print_motifs <- function(motif_list, n_letters){
  motifs <- unlist( lapply(motif_list, FUN = function(x){ print_single_motif(x, n_letters = n_letters) } ) )
  return(motifs)
}

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
