### For functions that still need a little polishing

base_features_from_signalsetlist <- function(x, section="interval", returns = "indices", wraptogRanges = "FALSE"){
  
  set_valleys <- positions_from_signalsetlist(x, "valleys", "indices")
  set_peaks <- positions_from_signalsetlist(x, "peaks", "indices")
  ## Needs optimizing, implement an apply function that can obtain
  ## All metadata in signalset that's not the base signal
  set_upstream <- lapply(x, '[[', 'distance_to_nearest_upstream_peak')
  set_downstream <- lapply(x, '[[', 'distance_to_nearest_downstream_peak')
  set_chr <-lapply(x, '[[', 'chromosome')
  set_start <-lapply(x, '[[', 'start')

  signalset_length <- length(x)
  base_feature_list <- vector(mode ="list",length = signalset_length)

  for(i in 1:signalset_length){
    ## Parse bw signal values associated to np interval
    signal <- x[[i]]$signal
    ## Access each observation
    peaks <- set_peaks[[i]]
    valleys <- set_valleys[[i]]
    
    ## Remove consecutive ocurrences
    peaks <- if(check_consecutive(peaks)==TRUE) remove_consecutive(peaks) else peaks
    valleys <- if(check_consecutive(valleys)==TRUE) remove_consecutive(valleys) else valleys
    
    sorted_index_vector <- sort(c(peaks,valleys))
    rep_sorted_index_vector <- c(sorted_index_vector[1],
                                 rep(sorted_index_vector[2:(length(sorted_index_vector)-1)],times=1,each=2),
                                 sorted_index_vector[length(unlist(sorted_index_vector))])
    
    base_features <- data.frame(start = rep_sorted_index_vector[seq(1,length(rep_sorted_index_vector),by=2)],
                              end = rep_sorted_index_vector[seq(2,length(rep_sorted_index_vector),by=2)])
    
    ## Calculate and Assign area base_features
    signal_area_length <- signal[rep_sorted_index_vector]
    base_features$chr <- set_chr[[i]]
    base_features$signal_area_rectangle_width <- diff(sorted_index_vector)
    base_features$signal_area_rectangle_height <- abs(rollapply(signal_area_length, 2, by=2, diff, partial = TRUE, align="left"))
    base_features$signal_rectangle_area <- (base_features$signal_area_rectangle_height * base_features$signal_area_rectangle_width)/2
    ## Add metadata
    ## Needs optimizing as well, same as above
    base_features$nearest_upstream <- set_upstream[[i]]
    base_features$nearest_downstream <- set_downstream[[i]]
    
    if(returns =="positions"){
      base_features$start <- base_features$start + (set_start[[i]]-1)
      base_features$end <- base_features$end + (set_start[[i]]-1)
      }

    if(section == "valley"){
      
      area_for_subset <- rep(base_features$signal_rectangle_area,each=2)
      width_for_subset <- rep(base_features$signal_area_rectangle_width, each=2)
      height_for_subset <- rep(base_features$signal_area_rectangle_height, each=2)
      matches <- match(rep_sorted_index_vector,valleys)
      ## Keep all valley values, discard peaks
      matches[is.na(matches)==FALSE] <- 1
      match_index <- matches
      match_index[!is.na(matches)] <- cumsum(match_index[!is.na(matches)])
      ## Assign 0, then subset as logical to keep valley values from area vector
      matches[is.na(matches)] <- 0
      ## Missing assertion to make valley vector and scores equal
      valley_area <- rollapply(area_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      valley_width <- rollapply(width_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      valley_height <- rollapply(height_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      
      if(returns == "positions") valleys <- valleys + (set_start[[i]]-1) else valleys
      
      base_feature_list[[i]] <- data.frame(valley = valleys, 
                                           chr = set_chr[[i]],
                                           area = valley_area,
                                           width = valley_width,
                                           height = valley_height,
                                           bps_to_next_np = rep(set_upstream[[i]], length(valleys)),
                                           bps_to_previos_np = rep(set_downstream[[i]], length(valleys)))
      }else{
        base_feature_list[[i]] <- base_features
      }
    
      
      
  }

  return(base_feature_list)

}
  



  
  
## Misc functions start here
# Useful for applying anything on rolling, fixed intervals. Imports Zoo package.
distances_between_points <- rollapply(sorted_index_vector, 2, by = 2, diff, partial = TRUE, align = "left")
signal_differences <- abs(rollapply(sorted_signals, 2, by = 2, diff, partial = TRUE, align = "left"))

## Loop used to determine if valleys or peaks should go first,
## Thanks to new implementation this is now deprecated.
if (valleys[[i]][1] > peaks[[i]][1]) {
  lead <- unlist(peaks)
  lagged <- unlist(valleys)
} else {
  lead <- unlist(valleys)
  lagged <- unlist(peaks)
}


## Join two equal sized vectors in an intertwining fashion
pos_index_vector <- c(rbind(lead,lagged))
pos_index_vector_i <- c(rbind(lead[diff(lead)>1],lagged))


##samplesig
x<-np_signals_from_bigwig(A549_chr1_bw,A549_ChIP_filtered[1:3])
