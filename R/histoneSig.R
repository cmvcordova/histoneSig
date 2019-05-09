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

check_and_load_refgenome <- function(refgenome){
  if(refgenome %in% BSgenome::available.genomes()){
    refgenome
  }else{
    stop("Provided reference genome not available.
         Consult available options with BSgenome::available.genomes()")}
  require(suppressPackageStartupMessages(refgenome),character.only = TRUE) # Load supplied genome
  genome <- eval(parse(text=refgenome)) ## Include in variable that will be used for further calls
  return(genome)
}

dna_one_hot <- function(x){
  
  # x is a basepair vector ex: 'A','C','G','T', etc
  ## Probably a horrible implementation, but it's doing the job
  ### HOW DO R FORMULA OBJECTS EVEN WORK?
  bases <- c('A','C','G','T')
  ## Build the encoder
  encoder = as.data.frame(bases)   
  names(encoder) <- 'bases'
  dmy <- caret::dummyVars(" ~ bases", data=encoder)
 ## Build one hot matrix from char vector
  x <- as.data.frame(x)
  names(x) <- 'bases'
  one_hot_matrix <- predict(dmy,x)
  colnames(one_hot_matrix) <- bases
  return(one_hot_matrix)
  
}

### SIGNALSET CLASS

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
  interpeak_right <- lead(start(np_object),1) - end(np_object) - 1
  interpeak_left <- c(NA,head(interpeak_right,-1))
  ## Remove NAs caused due to being the first or last peak in chromosome
  interpeak_right[is.na(interpeak_right)] <- 0
  interpeak_left[is.na(interpeak_left)] <- 0
  ## Set negative values caused by calculation at the start of each chromosome to 0
  interpeak_right[interpeak_right < 0] <- 0
  interpeak_left[interpeak_left < 0] <- 0
  ## Preallocate list to be populated with variable size signals
  signals <- vector(mode="list", length=n_peaks)
  signals <- lapply(1:n_peaks,function(x) signalSet())
  #########3 THINK ABOUT PREALLOCATING SIGNALSET'S SIGNAL WIDTH AS YOU ALREADY HAVE PRECOMPUTED WIDTHS
  ########## WHICH CORRESPOND TO THE SIGNAL'S LENGTH

    ## Fill rest of values in signalSet object
  for (i in (1:n_peaks)){

    signals[[i]]$signal <- granges_to_continuous(overlapper[lead_range_vector[i]:lagging_range_vector[i]])
    signals[[i]]$chromosome <- as.character(seqnames(np_object[i]))
    signals[[i]]$start <- starts[i]
    signals[[i]]$end <- ends[i]
    signals[[i]]$width <- widths[i]
    signals[[i]]$distance_to_nearest_upstream_peak <-   interpeak_right[i]
    signals[[i]]$distance_to_nearest_downstream_peak <- interpeak_left[i]

  }

  return(signals)

}


filter_signalSet <- function(signalsetlist, window_size, filter_function, ...){
  ## Use window_size for equal sized lowpass, fractional for proportional filter
  ## optional arguments for function
  ## TODO: add unused protection for dots
  ## BREAKS IF ALL SIGNALS IN THE SIGNALSET ARE THE SAME SIZE
  ## Doesn't work with singles
  ## Possibly fixed with the simplified = FALSE call of mapply (How to introduce it with variable function arguments present?)
  if(missing(...) == FALSE){
    dots <- list(...)
  for(i in 1:length(dots)) {
    assign(x = names(dots)[i], value = dots[[i]])
  }
    if('fractional' %in% names(dots)){
      window_size = ((unlist(lapply(signalsetlist,'[[', 'width'))) / fractional)
    }
  }
  ## Defaults to a lowpass filter
  if(missing(filter_function)) filter_function <- lowpass_filter else filter_function

  if(length(window_size) > 1){
  filtered_signals <- mapply(filter_function, x = lapply(signalsetlist,'[[', 'signal'), n = window_size)
  }else if(length(window_size) == 1){
  filtered_signals <- lapply(lapply(signalsetlist,'[[', 'signal'), filter_function, n = window_size)
  }
  else{stop("n must be a positive numeric integer or vector")}
  ## switch TS object back to numeric signal
  filtered_signals <- lapply(filtered_signals, as.numeric)
  # Replace NAs from timeseries with first real value in new, filtered vector
  ## find first non NA value in filter to subset later
  first_real_index <- lapply(filtered_signals,function(x) which(!is.na(x)))
  first_value_index <- lapply(first_real_index,min)
  values <- Map('[', filtered_signals, first_value_index)
  ##
  filtered_signals <- mapply(function(x,v) replace(x,is.na(x),v), x=filtered_signals ,v=values)
  for (i in 1:length(signalsetlist)){
  ## assign back into list
    ## Is there a way to do this with an apply or map?
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

  ## Obtained from a biostars post

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

base_features_from_signalsetlist <- function(x, section = "interval", returns = "indices",
                                             wraptoGRanges = "FALSE", unwrap = "FALSE",
                                             max_area_filter = "FALSE"){
  
  set_peaks <- positions_from_signalsetlist(x, "peaks", "indices")
  set_valleys <- positions_from_signalsetlist(x, "valleys", "indices")
  ## Rationale: even if there's n valleys, no area can be calculated without
  ## Accompanying peaks; therefore, these observations are dropped.
  no_peaks <- unlist(lapply(set_peaks,function(x) any(length(x)==0)))
  no_valleys <- unlist(lapply(set_valleys,function(x) any(length(x)==0)))
  ## Join into single logic vector for removal of observations with 0
  ## On peaks or valleys
  to_remove <- Reduce(`|`, list(no_peaks,no_valleys))
  ## subset valleys and peaks without 0 and continue
  set_valleys <- set_valleys[!to_remove]
  set_peaks <- set_peaks[!to_remove]
  x <- x[!to_remove]
  ## Needs optimizing, implement an apply function that can obtain
  ## All metadata in signalset that's not the base signal
  set_upstream <- lapply(x, '[[', 'distance_to_nearest_upstream_peak')
  set_downstream <- lapply(x, '[[', 'distance_to_nearest_downstream_peak')
  set_chr <-lapply(x, '[[', 'chromosome')
  set_start <-lapply(x, '[[', 'start')
  set_maximums <- sapply(lapply(x,'[[','signal'),max)

  signalset_length <- length(x)
  base_feature_list <- vector(mode ="list",length = signalset_length)

  for(i in 1:signalset_length){
    ## Parse bw signal values associated to np interval
    signal <- x[[i]]$signal
    ## Access each observation
    peaks <- set_peaks[[i]]
    valleys <- set_valleys[[i]]

    ## Remove consecutive ocurrences (is this still necessary?)
    peaks <- if(check_consecutive(peaks)==TRUE) remove_consecutive(peaks) else peaks
    valleys <- if(check_consecutive(valleys)==TRUE) remove_consecutive(valleys) else valleys

    sorted_index_vector <- sort(c(peaks,valleys))
    rep_sorted_index_vector <- c(sorted_index_vector[1],
                                 rep(sorted_index_vector[2:(length(sorted_index_vector)-1)],times=1,each=2),
                                 sorted_index_vector[length(unlist(sorted_index_vector))])
    
    start <- rep_sorted_index_vector[seq(1,length(rep_sorted_index_vector),by=2)]
    end <- rep_sorted_index_vector[seq(2,length(rep_sorted_index_vector),by=2)]
    extension <- diff(sorted_index_vector)

    base_features <- data.frame(start = start, end = end)

    ## Calculate and Assign area base_features
    signal_area_length <- signal[rep_sorted_index_vector]
    base_features$chr <- set_chr[[i]]
    base_features$extension <- extension
    base_features$height <- abs(zoo::rollapply(signal_area_length, 2, by=2, diff, partial = TRUE, align="left"))
    base_features$area <- (base_features$height * base_features$extension)/2
    
    
    ## Add metadata
    ## Needs optimizing as well, same as above
    base_features$signal_maximum <- set_maximums[i]
    base_features$bps_to_next_peak <- set_upstream[[i]]
    base_features$bps_to_previous_peak <- set_downstream[[i]]

    if(returns =="positions"){
      base_features$start <- base_features$start + (set_start[[i]]-1)
      base_features$end <- base_features$end + (set_start[[i]]-1)
    }

    if(section == "valley"){

      area_for_subset <- base_features$area
      extension_for_subset <- base_features$extension
      height_for_subset <- base_features$height
      matches <- match(rep_sorted_index_vector,valleys)
      groups <- matches
      ## Keep all valley values, discard peaks
      matches[is.na(matches)==FALSE] <- 1
      match_index <- matches
      match_index[!is.na(matches)] <- cumsum(match_index[!is.na(matches)])
      ## Assign 0, then subset as logical to keep valley values from area vector
      ## Missing assertion to make valley vector and scores equal
      ## calculate valley areas, extensions and heights on observations that contain
      ## a valley between two peaks
      valley_area <- as.vector(tapply(area_for_subset[match_index],groups,sum,na.rm=TRUE)) ## Area
      valley_extension <-  as.vector(tapply(extension_for_subset[match_index],groups,sum,na.rm=TRUE)) ## Extension
      valley_height <-  as.vector(tapply(height_for_subset[match_index],groups,sum,na.rm=TRUE)) ## Height

      if(returns == "positions") valleys <- valleys + (set_start[[i]]-1) else valleys

      base_feature_list[[i]] <- data.frame(valley = valleys,
                                           chr = rep(set_chr[[i]],length(valleys)),
                                           extension = valley_extension,
                                           height = valley_height,
                                           area = valley_area,
                                           signal_maximum = set_maximums[i],
                                           bps_to_next_peak = rep(set_upstream[[i]], length(valleys)),
                                           bps_to_previous_peak = rep(set_downstream[[i]], length(valleys)))

    } else {base_feature_list[[i]] <- base_features}
    ## Remove valleys with area = 0
    base_feature_list[[i]] <- base_feature_list[[i]][base_feature_list[[i]]$area!=0,]
    
    if(max_area_filter == "TRUE"){
      base_feature_list[[i]] <- base_feature_list[[i]][which.max(base_feature_list[[i]]$area),]
    }
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

  if(wraptoGRanges == "TRUE" & unwrap == "TRUE"){
  base_feature_list <- suppressWarnings(sort(do.call(c, unlist(base_feature_list,recursive=FALSE))))
  } else if(unwrap == "TRUE"){
  base_feature_list <- bind_rows(base_feature_list)
  } else {base_feature_list}

  return(base_feature_list)


}

build_feature_table <- function(x, metadata_as_features = FALSE, include_sequence = FALSE,  refgenome){
  ### To do: Generalize; merge this function with signal_matrix_from_signalSet to make it able to receive
  ## signalSets or Granges, whichever is supplied.
  
  ## This check must go, should be a feature of the experimental parsing function
  ### features from signalsetlist doesn't give return a signalset, therefore this check
  ## Returns an error on non granges base_feature_list
  
  if(class(x) != "GRanges" & class(x) != "list"){
    stop('Must provide a valid GRanges or signalSet list')
  } else if(class(x) == "list"){
    if(all(sapply(x, class) == 'signalSet') == "FALSE"){
      stop("List provided has at least one none signalSet object")}
  }
  
  if(include_sequence == TRUE){
    genome <- check_and_load_refgenome(refgenome)
  }
  
  ## Is there a better way to stop doing this at the beginning of each signalset associated function?
  ## A signalset accessor of sorts?
  
  if(class(x) == 'list'){
    
    starts <- sapply(signalSet,'[[','start')
    ends <- sapply(signalSet,'[[','end')
    chrs <- sapply(signalSet,'[[','chromosome')
    widths <- sapply(signalSet, '[[', 'width')
    chromosome <- unlist(mapply(rep,chrs,widths),use.names=FALSE)
    master_index <- unlist(mapply(`:`, starts, ends), use.names=FALSE)
    master_table <- data.table(master_index = master_index, chromosome = chromosome, signal = unlist(sapply(signalSet,'[[','signal')))
  } else if(class(x) == "GRanges"){
    ## Build an index to parse the aggregate sequence of all GRanges objects.
    master_grange <- x ## Copying x in case metadata has to be called later (which it does)
    elementMetadata(master_grange) <- NULL 
    master_granges_seqnames <- rep(as.vector(seqnames(master_grange)),width(master_grange))  ## obtain chrnames for index
    master_index <- unlist(mapply(`:`, start(master_grange), end(master_grange)), use.names=FALSE)   ## Create main bp vector
    # Store as a matrix to add further sections, depending on user supplied preferences
    master_table <- data.table(master_index, chromosome=master_granges_seqnames)}
  
  if(include_sequence == TRUE){
    ### to do: add option for custom genomes
    master_index_seqs <- as.character(getSeq(genome, master_grange))   ## Parse before deconstructing vector
    master_index_seqs <- unlist(strsplit(paste(master_index_seqs,collapse=""),"")) ## Split and join into per basepair sequence vector
    
    ## one hot encoding of previously obtained sequences
    master_seq_matrix <- dna_one_hot(master_index_seqs)
    master_table <- data.table(master_table,master_seq_matrix) ## Join into index + seq matrix
  }
  
  if(metadata_as_features == TRUE){
    metadata_features <- as.data.table(elementMetadata(x))
    metadata_features <- metadata_features[rep(seq_len(nrow(metadata_features)), width(x)),]
    master_table <- data.table(master_table,metadata_features)
  }
  
  return(master_table)
}


###### Experimental | Trying to build a parser to stop:
## 1) Is it a signalSet list or Signalset? discrepancy
## 2) stop the sapply line every single time we need an
## Attribute from our signalSet
### Todo: make to_parse able to receive multiple arguments.

parse_signalSet <- function(x, to_parse){

  ## Check if we're working with a full list of signalSet
  if(class(x) != "list" & class(x) != "signalSet"){
    stop('Must provide a valid signalSet list or signalSet to parse')
  } else if (all(sapply(x, class) == 'signalSet') == "FALSE"){
    stop("List provided has at least one none signalSet object")
  }

  force(to_parse)

  if(is.list(signalSet) == TRUE){
    parser <- function(to_parse) sapply(x,'[[', to_parse)
  } else
    parser <- function(to_parse) x$to_parse

  return(parser(to_parse))
}


per_chr_valley_plots <- function(x, separator='experiment', feature_to_plot = NULL ,
                                 condition=NULL, object_to_return="plot"){
  ## needs a bit more attention, currently designed
  ## to work with granges but a feature table approach should be in order
  ## To do: add argument to make something other than a segmented boxplot?
  if(class(x) != "GRanges" & class(x) != "data.frame"){
    stop('Must provide a valid GRanges, data.frame, or data.table')
  }
  
  if(!is.null(feature_to_plot)){
    if(!is.vector(feature_to_plot) & class(feature_to_plot) != "character"){
      stop('Must provid a valid character vector with the column names of the features that are to be plotted')
    } 
  }
  
  if(class(x) == "GRanges"){
    features <- as.data.table(elementMetadata(x))
    chrnames <- as.character(seqnames(x))
    ## 'positions' currently works since it's understood a valley object
    ## has the same start and end with length 1; this will
    ## undoubtedly break with a different sized input
    positions <- start(ranges(x)) 
    chr_data <- data.table(chrnames,positions,features) # separator subset features[,separator,with=FALSE]
  }else if(class(x) =="data.table"){
    chr_data <- data.table(x) ## for performance ?
  }
  ## To do: add a check between features_to_plot and
  ## names(data) to make sure something'll get plotted
  plot_names <- unique(chr_data$chrnames)
  n <- length(plot_names)
  ncol <- floor(sqrt(n))
  plot_list <- vector(mode = "list", length = n)
  
  for(i in 1:n){
    plotdata <- chr_data[chr_data[,chrnames == plot_names[i]]]
    
    plot_list[[i]] <- ggplot(data = plotdata, aes(x = chrnames, 
                                                  y = eval(parse(text = feature_to_plot)), 
                                                  fill = eval(parse(text = separator)))) +
      geom_boxplot() +
      xlab(plot_names[i]) +
      ylab(feature_to_plot) +
      guides(fill=guide_legend(title=separator)) +
      scale_y_log10() ## These transformations should be passed in as options to the function
    
  }
  names(plot_list) <- plot_names
  if(object_to_return == "plotlist"){
    return(plot_list)
  } else {
  return(do.call("grid.arrange", c(plot_list,ncol=ncol)))
  dev.off()
  }
  
}
