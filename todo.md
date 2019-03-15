Immediate  

- Get a robust, working implementation of signal_feature_matrix  
- Investigate why functions break when equal signals are provided (mapply), and
  when single observation signalsets are provided.
- Add parameter to signal_feature_matrix for different genomes based on BSgenome

Can wait a bit 

- Wrap signal_feature_matrix to genomic ranges  
- PONDER ABOUT adding an additional slot to the signalset class, "from
  experiment". This would have the goal of being able to feed a signalsetlist
to signal_feature_matrix and automatically create columns per experiment.  OR a
way to join only the experiment-exclusive signal column from
signal_feature_matrix objects (Figure out which is more or less computationally
expensive).

Long Term

- Add functionality to plotSignal to highlight regions that overlap regions of
  interest (ex: If plotting ChIP signal, show which ranges overlap with ATAC if
predicting open chromatin).
- When positions instead of indices are specified in plotSignal, adjust x axis
  to show this.
- See if theres an object that can compactly store signals 
- If above fails, HDF5 File compatibility for dampening the cost associated to
  representing signals in a long format
- Build a wrapper for being able to read multiple file pairs for cell/tissue
  type experimenting.

Thoughts that could probably culminate in some useful implementation:

- Can you build a general accesor function or redesign signalset classes to
  stop having to lapply every single parameter at the start of a function?
