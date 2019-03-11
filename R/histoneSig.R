# Place in required: rtracklayer, GenomicRanges, gtools, dplyr, tidyverse, quantmod
### Include function for creating Homo.sapiens.hg38 database from UCSC track

## HELPER FUNCTIONS

import.np <- function(np_file){
  ## Declare columns to accomodate narrowPeak format
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  ## Read in calling import with extra parameters
  return(import(np_file,format="bed",extraCols=extraCols_narrowPeak))
}

## base filtering function; imported from github (reference pending)
lowpass_filter <- function(x,n){stats::filter(x,rep(1/n,n), sides=1)}

check_consecutive <- function(x){
  any(rle(diff(x))$values == 1)
}

remove_consecutive <- function(x){
  ## Should look for an in-place solution
  split_by_consecutive_intervals <- split(x, cumsum(c(1, diff(x) != 1)))
  consecutive_intervals <- split_by_consecutive_intervals[rapply(split_by_consecutive_intervals,length,how="list") > 1]
  consecutive_intervals <- lapply(consecutive_intervals,function(x) min(x) + (max(x) - min(x) +1)/2)
  return(unlist(modifyList(split_by_consecutive_intervals,consecutive_intervals),use.names=FALSE))

}

## Function to filter and order genomic ranges by 23 canonical chromosomes
granges_chr_filter <- function(granges_object){
  granges_object <- granges_object[grep("[chr]+([0-9]{1,2}$|X|Y)",seqnames(granges_object)),]
  seqlevels(granges_object) <- unique(as.character(seqnames(granges_object)))
  granges_object <- sortSeqlevels(granges_object)
  granges_object <- sort(granges_object)
}


####### Same as function from quantmod, only difference
####### Is the sum is +1 instead of +2

findvalleysone <- function (x, thresh = 0)
{
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) >
                 0) + 1
  if (!missing(thresh)) {
    if (sign(thresh) > 0)
      thresh <- -thresh
    pks[x[pks - 1] - coredata(x[pks]) < thresh]
  }
  else pks
}

findpeaksone <- function (x, thresh = 0)
{
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) <
                 0) + 1
  if (!missing(thresh)) {
    if (sign(thresh) < 0)
      thresh <- -thresh
    pks[x[pks - 1] - coredata(x[pks]) > thresh]
  }
  else pks
}

#Obtain position from ranges

granges_to_continuous <- function(x){
  #Where x is a GenomicRanges object with an associated signal value
  #(parsed from bigwig)
  score_vector <- rep(x$score,width(x))
  return(score_vector)
}

signalSet <- function(){

  signal <- list(signal = vector(),
                 chromosome = character(),
                 start = numeric(),
                 end = numeric(),
                 width = numeric() ,
                 distance_to_nearest_upstream_peak = numeric(),
                 distance_to_nearest_downstream_peak = numeric())

  class(signal) <- "signalSet"

  signal
}

signals_from_bigwig <- function(bw_object){
######### currently only as a test function, use only at own risk
  ##Make object that will be transformed to signal given a bigwig file

  widths <- width(bw_object)
  starts <- start(bw_object)
  ends<- end(bw_object)
  n <- length(bw_object)
  ## calculate interpeak distance, set right and left orientation
  ### THIS MUST BE FIXED PER CHROMOSOME
  interpeak_right <- lead(starts,1) - ends - 1
  interpeak_left <- c(NA,head(interpeak_right,-1))
  ## Preallocate list to be populated with variable size signals
  signals <- vector(mode="list", length=n)
  signals <- lapply(1:n,function(x) signalSet())
  #########3 THINK ABOUT PREALLOCATING SIGNALSET'S SIGNAL WIDTH AS YOU ALREADY HAVE PRECOMPUTED WIDTHS
  ########## WHICH CORRESPOND TO THE SIGNAL'S LENGTH
  ## maybe https://www.r-bloggers.com/how-to-use-lists-in-r/ for handling lists

  ## Fill rest of values in signalSet object
  for (i in (1:n)){

    signals[[i]]$signal <- granges_to_continuous(bw_object[i])
    signals[[i]]$chromosome <- as.character(seqnames(bw_object[i]))
    signals[[i]]$start <- starts[i]
    signals[[i]]$end <- ends[i]
    signals[[i]]$width <- widths[i]
    signals[[i]]$distance_to_nearest_upstream_peak <-   interpeak_right[i]
    signals[[i]]$distance_to_nearest_downstream_peak <- interpeak_left[i]

  }

  return(signals)

}

np_signals_from_bigwig <- function(bw_object, np_object){

  ##Make object that will be transformed to signal given specified peakfile, bigwig pairs

  ## Calculate overlap vectorS
  overlapper <- subsetByOverlaps(bw_object,np_object)
  overlaps <- countOverlaps(np_object, bw_object)
  widths <- width(np_object)
  starts <- start(np_object)
  ends<- end(np_object)
  ## create lead and lagging vectors to subset overlap object
  n_peaks <- length(overlaps)
  lagging_range_vector <- cumsum(overlaps)
  lead_range_vector<- c(1,lagging_range_vector+1)[1:n_peaks]
  ## calculate interpeak distance, set right and left orientation
  ### THIS MUST BE FIXED PER CHROMOSOME
  interpeak_right <- lead(start(np_object),1) - end(np_object) - 1
  interpeak_left <- c(NA,head(interpeak_right,-1))
  ## Preallocate list to be populated with variable size signals
  signals <- vector(mode="list", length=n_peaks)
  signals <- lapply(1:n_peaks,function(x) signalSet())
  #########3 THINK ABOUT PREALLOCATING SIGNALSET'S SIGNAL WIDTH AS YOU ALREADY HAVE PRECOMPUTED WIDTHS
  ########## WHICH CORRESPOND TO THE SIGNAL'S LENGTH
  ## maybe https://www.r-bloggers.com/how-to-use-lists-in-r/ for handling lists

  ## Fill rest of values in signalSet object
  for (i in (1:n_peaks)){

    signals[[i]]$signal <- granges_to_continuous(overlapper[lead_range_vector[i]:lagging_range_vector[i]])
    signals[[i]]$chromosome <- as.character(seqnames(overlapper[i]))
    signals[[i]]$start <- starts[i]
    signals[[i]]$end <- ends[i]
    signals[[i]]$width <- widths[i]
    signals[[i]]$distance_to_nearest_upstream_peak <-   interpeak_right[i]
    signals[[i]]$distance_to_nearest_downstream_peak <- interpeak_left[i]

  }

  return(signals)

}


lowpass_filter_signalSet <- function(signalsetlist,n){

  ## There must be a better way to do this
  ## workaround to R's filter function addition of a time series message
  ### See if this is removable within the package

  filtered_signals <- lapply(lapply(signalsetlist,'[[', 1), function(x){lowpass_filter(x,n)})
  filtered_signals <- lapply(filtered_signals, as.numeric)

  for (i in 1:length(signalsetlist)){

    signalsetlist[[i]]$signal <- filtered_signals[[i]]
  }

  return(signalsetlist)
}

## Where query is the file to obtain regions with multiple overlaps from
## and target is the file that can overlap multiple times with query

parse_regions_without_multiple_overlaps <- function(query, target){
  overlap_vector <- countOverlaps(query,target)
  query$num_overlapping_ranges <- overlap_vector
  query[countOverlaps(query,target) <= 1,]
}


## This is used for informative purposes only, so as to subset negative examples
parse_regions_with_multiple_overlaps <- function(query, target){
  overlap_vector <- countOverlaps(query,target)
  query$num_overlapping_ranges <- overlap_vector
  query[countOverlaps(query,target) > 1,]
}


############ TEST FURTHER FOR ROBUSTNESS
############ AND GENERAL SENSE
### Probably rename into "join overlaps" or something

obtain_extended_parsing_regions <- function(query,target){

  #Given a query and a target, returns an object with target values
  #that joins regions that share overlap regions in the corresponding query

  ## Subset target's zeros out
  target$num_overlaps <- countOverlaps(target,query)
  target_length <- length(target)
  zeros <- target[target$num_overlaps == 0]
  #Obtain ranges with at least one overlap
  target_overlap_ranges <- target[target$num_overlaps > 0]

  which <- findOverlaps(target_overlap_ranges,query)
  ## Factor
  factored_indices <- factor(to(which))

  # parsing_index allows us to see which were the query's overlapping ranges with the target
  parsing_index <- from(which)
  target_ranges <- target_overlap_ranges[parsing_index,]
  ## an if will have to be added here to prevent addition of ranges if they're
  ## longer than a threshold
  iterator <- data.frame(index = unlist(lapply(split(factored_indices,factored_indices),seq_along)),
                         target_start = start(target_ranges),
                         target_end = end(target_ranges),
                         width = width(target_ranges),
                         query_overlap_range = factored_indices,
                         parsing_index = parsing_index,
                         chr = seqnames(target_ranges))

  ## probably not the most efficient way
  pre_ranges <- split(iterator, f = iterator$query_overlap_range)
  new_start <-  sapply(pre_ranges,function(x) x[1,2])
  new_end <- sapply(pre_ranges,function(x) x[length(x[,3]),3])
  chr <- sapply(pre_ranges,function(x) x[1,7])

  ## change chr value to dynamically adjust according to chr experiment that's provided
  new_multiples <- GRanges(seqnames = chr, ranges = IRanges(start = new_start,end = new_end))
  elementMetadata(zeros) <- NULL
  new_object <- c(zeros,new_multiples)
  new_object$intervalwidth <- width(new_object)
  new_object$num_overlaps <- countOverlaps(new_object,query)

  return(unique(sort(new_object)))
  ## further optimization: remove ranges shared between
  ## created range object and within-range-object-singles (if any)

}

grextend <- function(x, upstream=0, downstream=0){

  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

positions_from_signalsetlist <- function(signalsetlist, points, returns = "positions"){
  ##### Breaks on non-nested lists (working as intended?)
  ##### Search for default methods on S3 classes
  #### create methods for both signalSet and signalSetlist
  ##(or don't, learn to circumvent the list wrapper in signalSets
  force(points)
  if(points == "valleys"){
    x <- lapply(lapply(signalsetlist,'[[', 1), function(x){findvalleysone(x)})
  }else if(points == "peaks"){
    x <- lapply(lapply(signalsetlist,'[[', 1), function(x){findpeaksone(x)})
  }else{
    stop("points argument must be specified as either peaks or valleys")
  }
  ##apply findValleys to return alleys in sets
  ## Obtain respective starting positions, subtract 1 since these will be
  ## added to the valley's index
  if (returns == "indices"){
    return(x)
  }else{
    ### subtracting 1 from starts to report valleys correctly after sum
    starts <- unlist(lapply(signalsetlist,'[[', 3))-1
    return(Map(`+`, x, starts))
  }

}


## Base plotting function
plotSignal <- function(x, highlight="none", ...){
  ### Receives a single signalSet
  ### plotting function to aid in the illustration of signalsets
  ### will show maxima and/or minima as specified with lines crossing each index point
  ### ideally, it will also show areas.
  signal = x[[1]]$signal
  ## add, if NOT list just take as is

  dt = data.table(position = seq_along(signal), signal = signal)
  setkey(dt,position)

  p = ggplot(data = dt , aes(x = position, y= signal)) +
    labs(x="Position",y="Signal Value") +
    geom_line()

  if(highlight == "valleys"){

    valleys <- data.table(position = unlist(positions_from_signalsetlist(x, "valleys", "indices")))
    dt[valleys, on = "position", valley_exists := i.position]
    p = p + geom_point(data= dt[dt$valley_exists>0], color="red") +
      scale_color_discrete(name = "Point", labels = "Valleys")

  } else if(highlight == "peaks") {

    peaks <- data.table(position = unlist(positions_from_signalsetlist(x, "peaks", "indices")))
    dt[peaks, on = "position", peak_exists := i.position]
    p = p + geom_point(data= dt[dt$peak_exists>0], color="blue")

  } else if(highlight =='both') {

    valleys <- data.table(position = unlist(positions_from_signalsetlist(x, "valleys", "indices")))
    peaks <- data.table(position = unlist(positions_from_signalsetlist(x, "peaks", "indices")))
    dt[valleys, on = "position", valley_exists := i.position]
    dt[peaks, on = "position", peak_exists := i.position]
    p = p + geom_point(data= dt[dt$valley_exists>0], color="red") +
      geom_point(data= dt[dt$peak_exists>0], color="blue")

  }

  return(p)
}

signal_area_descriptors_from_signalsetlist <- function(x, section="interval"){

## Not working with "positions" argument instead of indices
  set_valleys <- positions_from_signalsetlist(x, "valleys", "indices")
  set_peaks <- positions_from_signalsetlist(x, "peaks", "indices")

  signalset_length <- length(x)
  descriptor_list <- vector(mode ="list", length = signalset_length)

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

    descriptors <- data.frame(start = rep_sorted_index_vector[seq(1,length(rep_sorted_index_vector),by=2)],
                              end = rep_sorted_index_vector[seq(2,length(rep_sorted_index_vector),by=2)])

    signal_area_length <- signal[rep_sorted_index_vector]
    descriptors$signal_area_rectangle_width <- diff(sorted_index_vector)
    descriptors$signal_area_rectangle_length <- abs(rollapply(signal_area_length, 2, by=2, diff, partial = TRUE, align="left"))
    descriptors$signal_rectangle_area <- (descriptors$signal_area_rectangle_length * descriptors$signal_area_rectangle_width)/2

    if(section == "valley area"){

      area_for_subset <- rep(descriptors$signal_rectangle_area,each=2)
      matches <- match(rep_sorted_index_vector,valleys)
      ## Keep all valley values, discard peaks
      matches[is.na(matches)==FALSE] <- 1
      match_index <- matches
      match_index[!is.na(matches)] <- cumsum(match_index[!is.na(matches)])
      ## Assign 0, then subset as logical to keep valley values from area vector
      matches[is.na(matches)] <- 0
      ## Missing assertion to make valley vector and scores equal
      valley_scores <- rollapply(area_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      descriptor_list[[i]] <- data.frame(valley = valleys, value = valley_scores)

    }else if(section == "interval"){

      descriptor_list[[i]] <- descriptors

    }
  }
  if(returns == "indices"){
    return(descriptor_list)
  }else if(returns == "positions"){


}

}

## Function useful for determining baseline percentage of overlaps
## Between a given query and a target
overlap_baseline <- function(query, target, return_unique = "FALSE"){
  ## Warning: has no protection for cases where query and target have different
  ## Seqnames; use function with this taken into consideration.
  ## Currently, user must remove unwanted chromosomes from granges object and make
  ## Sure they correspond in both datasets

  overlaps = countOverlaps(query,target)
  if(return_unique == "TRUE") overlaps[overlaps!=0] <- 1 else overlaps

  query_ranges_by_chr <- tibble(chr = as.vector(seqnames(query)), overlaps = overlaps) %>%
    group_by(chr) %>%
    summarize_all(sum) %>%
    add_row(., chr = "total", overlaps = sum(.$overlaps))


  target_ranges_by_chr <- tibble(chr = as.vector(seqnames(target)), count = rep(1,length(target))) %>%
    group_by(chr) %>%
    summarize_all(sum) %>%
    add_row(., chr = "total", count = sum(.$count))

  overlap_tibble <- tibble(chr = target_ranges_by_chr$chr,
                           overlap_percentage = (unlist(query_ranges_by_chr$overlaps) /
                                                   unlist(target_ranges_by_chr$count)) *100)
  return(overlap_tibble)

}


## Still prototyping with it, use at own risk
## Width calculation is off.

base_features_from_signalsetlist <- function(x, section="interval", returns = "indices", wraptoGRanges = "FALSE"){

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
    base_features$extension <- diff(sorted_index_vector)
    base_features$height <- abs(rollapply(signal_area_length, 2, by=2, diff, partial = TRUE, align="left"))
    base_features$area <- (base_features$height * base_features$extension)/2
    ## Add metadata
    ## Needs optimizing as well, same as above
    base_features$bps_to_next_peak <- set_upstream[[i]]
    base_features$bps_to_previous_peak <- set_downstream[[i]]

    if(returns =="positions"){
      base_features$start <- base_features$start + (set_start[[i]]-1)
      base_features$end <- base_features$end + (set_start[[i]]-1)
    }

    if(section == "valley"){

      area_for_subset <- rep(base_features$area,each=2)
      extension_for_subset <- rep(base_features$extension, each=2)
      height_for_subset <- rep(base_features$height, each=2)
      matches <- match(rep_sorted_index_vector,valleys)
      ## Keep all valley values, discard peaks
      matches[is.na(matches)==FALSE] <- 1
      match_index <- matches
      match_index[!is.na(matches)] <- cumsum(match_index[!is.na(matches)])
      ## Assign 0, then subset as logical to keep valley values from area vector
      matches[is.na(matches)] <- 0
      ## Missing assertion to make valley vector and scores equal
      valley_area <- rollapply(area_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      valley_extension <- rollapply(extension_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")
      valley_height <- rollapply(height_for_subset[as.logical(matches)], 2, by = 2, sum, partial = TRUE, align = "left")

      if(returns == "positions") valleys <- valleys + (set_start[[i]]-1) else valleys

      base_feature_list[[i]] <- data.frame(valley = valleys,
                                           chr = set_chr[[i]],
                                           extension = valley_extension,
                                           height = valley_height,
                                           area = valley_area,
                                           bps_to_next_peak = rep(set_upstream[[i]], length(valleys)),
                                           bps_to_previous_peak = rep(set_downstream[[i]], length(valleys)))

    } else {base_feature_list[[i]] <- base_features}
  }

  if(wraptoGRanges == "TRUE"){
    ## Simplify valley dataframe names before wrapping
    ## as GRanges object
    if(section == "valley"){
      base_feature_list <- lapply(base_feature_list, function(x){
        names(x)[names(x) == "valley"] <- "start"
        x$end <- x$start
        return(x)
      })
    }

    base_feature_list <- lapply(base_feature_list, function(x){
      makeGRangesFromDataFrame(x, seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end",
                               keep.extra.columns=TRUE)
    })

  }

  return(base_feature_list)


}
