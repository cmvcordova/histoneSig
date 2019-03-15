Immediate  

- Get a robust, working implementation of signal_feature_matrix  


Can wait a bit 

- Wrap signal_feature_matrix to genomic ranges  
- PONDER ABOUT adding an additional slot to the signalset class, "from
  experiment". This would have the goal of being able to feed a signalsetlist to signal_feature_matrix and automatically create columns per experiment. 

Long Term

- Add functionality to plotSignal to highlight regions that overlap regions of
  interest (ex: If plotting ChIP signal, show which ranges overlap with ATAC if
predicting open chromatin).
- See if theres an object that can compactly store signals 
- If above fails, HDF5 File compatibility for dampening the cost associated to
  representing signals in a long format
