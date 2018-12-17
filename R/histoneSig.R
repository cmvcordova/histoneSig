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

## Function to filter and order genomic ranges by 23 canonical chromosomes
granges_chr_filter <- function(granges_object){
  granges_object <- granges_object[grep("[chr]+([0-9]{1,2}$|X|Y)",seqnames(granges_object)),]
  seqlevels(granges_object,force=TRUE) <- unique(as.character(seqnames(granges_object)))
  granges_object <- sortSeqlevels(granges_object)
  granges_object <- sort(granges_object)
}

#Obtain position from ranges

granges_to_continuous <- function(x){
  #Where x is a GenomicRanges object with an associated signal value
  #(parsed from bigwig)
  score_vector <- rep(x$score,width(x))
  return(score_vector)
}

################################# CURRENTLY WORKING ON THIS ##################################

#### IN DEVELOPMENT

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


signals_from_bigwig <- function(bw_object, np_object){

  ##Make object that will be transformed to signal given specified peakfile, bigwig pairs

  ## Calculate overlap vector
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
    new_object$width <- width(new_object)
    new_object$num_overlaps <- countOverlaps(new_object,query)


    #### Create subobject to join later

    ##### Sum FROM subjecthit together on subjecthit TO vector

    ## further optimization: remove ranges shared between
    ## created range object and within-range-object-singles (if any)

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


### Subsetting with only numbers for S3 class is pretty weak, find alternate form
valley_positions_from_signalsetlist <-function(signalsetlist){
  ####### REALLY WEAK FUNCTION, REWRITE ASAP
  #### BREAKS WHEN A SINGLE SIGNALSET IS GIVEN
  ####### BREAKS ON NON NESTED LISTS
  ##### Search for default methods on S3 classes
  #### create methods for both signalSet and signalSetlist
  ##(or don't, learn to circumvent the list wrapper in signalSets

  ##apply findValleys to returnv alleys in sets
  valleys <- lapply(lapply(signalsetlist,'[[', 1), function(x){findValleys(x)})
  ## Obtain respective starting positions, subtract 1 since these will be
  ## added to the valley's index
  starts <- unlist(lapply(signalsetlist,'[[', 3))-1

  ### subtracting 1 from starts to report valleys correctly
  ### this should probably not happen here
  return(Map(`+`, valleys, starts))

}

## REMEMBER TO ADD findpeakspeaks fixed function (+1 instead of +2)

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

############ TEST FURTHER FOR ROBUSTNESS
############ AND GENERAL SENSE

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
  new_object$width <- width(new_object)
  new_object$num_overlaps <- countOverlaps(new_object,query)


  ## further optimization: remove ranges shared between
  ## created range object and within-range-object-singles (if any)

}

