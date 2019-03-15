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

Consult the sessionInfo.txt file within this repository to see which files are
being used in a development environment.  

```R
devtools::install_github('semibah/histoneSig')
```

## Demo of current utilities

importing bw,np
reading a signal + plot
signal filtering + plot
geometric feature extraction (basefeatures; show df +granges)
easy to use implementation for comparing between regions of generated signals from (basefeaturedf, showdf +granges)
different -seq experiments, against a reference genome or both. 

## Built With

* [RStudio](https://www.rstudio.com/) - Both desktop and server versions.
* [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) - For documentation purposes.

## Contributing

Nothing formal here, [just drop me a line](mailto:cesarmiguelv@gmail.com)

## Authors

* **Cesar Miguel Valdez** - *Initial work* - [semibah](https://github.com/semibah)

See also the list of [contributors](https://github.com/semibah/histoneSig/contributors) who participated in this project. Currently empty; you could be the first!

## License

Pending - probably an MIT one in time.

## Acknowledgments

* Ensenada Center for Scientific Research and Higher Education (CICESE) - Dr.
  Carlos Brizuela & Dr. Ivetth Corona, for their attentive guidance and
valuable suggestions.  
* National Council on Science and Technology (CONACyT) - Generously provided a
  scholarship.  
* Mexican Community of Bioinformatic Software Developers (CDSB) - Dr. Leonardo
  Collado-Torres & Dr. Alejandro Reyes, for welcoming me to CDSB, inspiring me
to create an Rpackage alongside my thesis and general guidance.  
* The internet - it can be a pretty nice place, after all.  

