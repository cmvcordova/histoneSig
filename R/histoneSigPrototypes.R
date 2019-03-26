### For functions that still need a little polishing


feature_matrix_from_signalSet <- function(signalSet){

}

# Parse granges objects as provided and join to main index.
#### make sub data tables and then join to main matrix
for(for i in x){}

master_seq_matrix[miniseq, V2 := i.V2, on =c("master_granges_index", "master_granges_seqnames")]
## miniseq is current sub matrix, V2 is a placeholder name for signal value col
rm(test_seq_matrix)
test_matrix <- na.omit(test_matrix)
setkey(test_matrix,V1)


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

## Subset observations that have "valley" in names
## Setting up start and end vectors to make GRange wrapping easier
lapply(base_feature_list[unlist(lapply(base_feature_list,
                                       function(x) ("valley" %in% names(x))==TRUE))],
       function(x){
         names(x)[names(x) == "valley"] <- "start"
         x$end <- x$start
       })


##samplesig
x<-np_signals_from_bigwig(A549_chr1_bw,A549_ChIP_filtered)
chr1_chips <- A549_ChIP_filtered[seqnames(A549_ChIP_filtered) == 'chr1']

x<-np_signals_from_bigwig(A549_chr1_bw,chr1_chips[1128:1130])

base_features_from_signalsetlist(z, section="valley", returns="positions", wraptoGRanges = "TRUE", unwrapGRanges = "TRUE")

