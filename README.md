# histoneSig

Currently a proof-of-concept tool. With the introduction of the signalSet class
and related methods, it allows for the generation and manipulation of
continuous, per-basepair, signal-like representations from -seq related bed and
bigwig files.   

Current utilities include signal filtering, basic geometric feature extraction
and an easy to use implementation for comparing between regions of generated
signals from different -seq experiments, against a reference genome or both.  

Features calculated from signals can then be integrated alongside
readily-provided -seq average values. Emphasizing the use of diverse sources of
information, our package interfaces with BSgenome for the integration of
sequence information within regions deemed significant by signal-obtained
metrics.  

While still in a prototype stage, histoneSig will ideally allow for downstream
analyses augmented by features obtained directly from our proposed signal
representation, wrapped in efficient data objects that will trivialize
interoperation with related analysis tools.  

## Installation

Get R 3.5.1.  

Consult the sessionInfo.txt file within this repository to see what's
being used in a development environment.  

```R

devtools::install_github('semibah/histoneSig')

```

## Demo of current utilities

Starting from a narrow peak (or any other bedfile) and its corresponding bigwig file, we're able to:

### Creating continuous signal objects from bedfile, bigwig pairs

```R 

## Load our peak or preferred bedfile np_file <-
import.np('path/to/npfile.bed')

## Set the ranges we just got obtained to parse relevant bigWig fragments
parsing_bw_ranges <- granges_chr_filter(np_file)

## Parse bigWig bw_file <- import.bw(con = BigWigFile("path/to/bwfile.bigWig"),
selection = BigWigSelection(parsing_bw_ranges))

## Obtain signals from both of our files 
your_first_signalset <- np_signals_from_bigwig(np_file, bw_file) 

```  

Behold, a `signalSet` observation in all its splendor  

![signalSet](/readmeimgs/signalsetclassexample.png)  

### Filter obtained signals as a set  
 
The default method is a lowpass filter. Said
filter can take a fixed `window_size` or an equal fraction of each signal in the
set as a `fractional` window. You can also pass your own filter functions to
`filter_signalSet()` (results may vary).

```R 

filtered_signalset <- filter_signalSet(your_first_signalset, fractional = 25) 

```

Now, let's use `plotSignal()` to compare the first signal of the `signalSet` we've
obtained.


```R

rawsignalplot <- plotSignal(your_first_signalset[1]) 
filteredsignalplot <- plotSignal(filtered_signalset[1])
gridExtra::grid.arrange(rawsignalplot, filteredsignalplot, ncol=2)

``` 
![signal plots](/readmeimgs/sidebysideplots.png)

We may also illustrate detected peaks (blue) and valleys (red). These will then
be used as references to calculate geometric features.

```R 

plotSignal(filtered_signalset[1], highlight="both") 

```
![filtered plot with highlights](/readmeimgs/filteredplothighlights.png)

### Calculate geometric features based on per-signal detected valleys and peaks

Calculating base features from a given `signalSet` is now possible; if posterior
interaction with `GenomicRanges` objects is desired, we can set our `wraptoGRanges`
argument as `TRUE`; else, we'll obtain a `data.frame`. Here, we'll specify notable
valleys found in our signal and their associated geometric features: valley
width ("extension"), height, area and distances to next and previous peaks in
the provided bedfile.

```R 

base_features_from_signalsetlist(filtered_signalset,
				 section="valley", returns="positions", wraptoGRanges=TRUE) 
``` 
![valley features as GRanges](/readmeimgs/basefeatures.png)

### Represent obtained features alongside signal values and reference sequences per genomic position

Finally, for comparative analyses, we may create a feature table from a
`GRanges` or `signalSet`. Sequence information may be integrated, as a one-hot matrix, setting
the `include_sequence` parameter to `TRUE`. We'll obtain a neat representation
which may or may not include neat signal and geometric features in its
metadata, as a `data.table`. This can then be easily interfaced with other
R libraries/models/packages. 

So, from a given `GRanges` (or vanilla `signalSet`) of the following kind:  

![generic GRanges object](/readmeimgs/genericgranges.png)

After running the following command,
```R 

build_feature_table(generic_GRanges, metadata_as_features = TRUE, include_sequence =
TRUE,  refgenome = "BSgenome.Hsapiens.UCSC.hg38")
 

```
We'd obtain a `data.table` like this

![feature table](/readmeimgs/genericfeaturetable.png)

## Built With

* [RStudio](https://www.rstudio.com/) - Both desktop and server versions.
* [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html)- For documentation purposes.

## Contributing

Nothing formal here just yet, [just drop me a line](mailto:cesarmiguelv@gmail.com)

## Authors

* **César Miguel Valdez Córdova** - *Initial work* - [semibah](https://github.com/semibah)

See also the list of
[contributors](https://github.com/semibah/histoneSig/contributors) who
participated in this project. Currently empty; you could be the first one!

## License

Pending - probably an MIT one in time.

## Acknowledgments

* Ensenada Center for Scientific Research and Higher Education (CICESE) - Dr.
  Carlos Brizuela & Dr. Ivetth Corona, for their attentive guidance and
valuable suggestions.  
* Mexican National Council on Science and Technology (CONACyT) - Generously
  provided a scholarship for the development of my master's thesis, which gave
way to histoneSig.  
* Mexican Community of Bioinformatic Software Developers (CDSB) - Dr. Leonardo
  Collado-Torres & Dr. Alejandro Reyes, for welcoming me to the CDSB, inspiring
me to create an R package alongside my thesis and general guidance.  
* The internet - it can be a pretty nice place, after all.  

