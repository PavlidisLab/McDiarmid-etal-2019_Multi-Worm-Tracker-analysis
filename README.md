# McDiarmid-etal-2019_Multi-Worm-Tracker-analysis
Analysis code for "Systematic phenomics analysis of ASD-associated genes reveals shared functions and parallel networks underlying reversible impairments in habituation learning" (DOI: 10.1101/687194)

<!-- TABLE OF CONTENTS -->
## Table of Contents
* [About the Project](#about-the-project)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
  * [Data](#data)
  * [Additional data](#additional-data)
* [Usage](#usage)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project

* Preprint (BioRxiv): [https://doi.org/10.1101/687194](https://doi.org/10.1101/687194)
* Data repository (Dataverse): [https://doi.org/10.5683/SP2/FJWIL8](https://doi.org/10.5683/SP2/FJWIL8)
* Analysis source code repository (Github): [https://github.com/PavlidisLab/McDiarmid-etal-2019_Multi-Worm-Tracker-analysis/](https://github.com/PavlidisLab/McDiarmid-etal-2019_Multi-Worm-Tracker-analysis/)

>SUMMARY:  
>A major challenge facing the genetics of Autism Spectrum Disorders (ASD) is the large and growing number of candidate risk genes and gene variants of unknown functional significance. Here, we used Caenorhabditis elegans to systematically functionally characterize ASD-associated genes in vivo. Using our custom machine vision system we quantified 26 phenotypes spanning morphology, locomotion, tactile sensitivity, and habituation learning in 87 strains each carrying a mutation in an ortholog of an ASD-associated gene. We identified hundreds of novel genotype-phenotype relationships ranging from severe developmental delays and uncoordinated movement to subtle deficits in sensory and learning behaviors. We clustered genes by similarity in phenomic profiles and used epistasis analysis to discover parallel networks centered on CHD8•chd-7 and NLGN3•nlg-1 that underlie mechanosensory hyper-responsivity and impaired habituation learning. We then leveraged our data for in vivo functional assays to gauge missense variant effect. Expression of wild-type NLG-1 in nlg-1 mutant C. elegans rescued their sensory and learning impairments. Testing the rescuing ability of all conserved ASD-associated neuroligin variants revealed varied partial loss-of-function despite proper subcellular localization. Finally, we used CRISPR-Cas9 auxin inducible degradation to determine that phenotypic abnormalities caused by developmental loss of NLG-1 can be reversed by adult expression. This work charts the phenotypic landscape of ASD-associated genes, offers novel in vivo variant functional assays, and potential therapeutic targets for ASD.

Please see the pre-print for more information about the project.


<!-- GETTING STARTED -->
## Getting Started

Make sure that you have R installed with the requirements listed below. Raw and processed data can be obtained from the Dataverse link mentionned above.

### Prerequisites
The code was executed using the following R version and platform:
```
R version 3.5.0 (2018-04-23) -- "Joy in Playing"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)
```
#### System packages
Some system packages are required for specific scripts. The operating system used to develop those scripts was `CentOS Linux release 7.6.1810 (Core)`. The following packages are required for specific parts of the document:

* The system package for Cairo Graphics is required to save some of the figures in pdf format. See instructions at [https://cairographics.org/download/](https://cairographics.org/download/).
* ImageMagick is used to process some images through system calls. The version tested was `ImageMagick 6.7.8-9 2016-06-16 Q16`. See http://www.imagemagick.org 

Please report any missing system packages by [opening an issue](https://github.com/PavlidisLab/McDiarmid-etal-2019_Multi-Worm-Tracker-analysis/issues) on this repository and mentionning which package was missing as well as which operating system was used.

### Installation
The R packages required for this project are listed in `requirements.R`. There may be additional system packages (i.e. Cairo Graphics) that need to be installed outside of the R environment.

Please report any missing dependencies/requirements by [opening an issue](https://github.com/PavlidisLab/McDiarmid-etal-2019_Multi-Worm-Tracker-analysis/issues) on this repository.

### Data
Raw or preprocessed data can be obtained from our lab dataverse. See: `datasets/README.md` or `processed_data/README.md`.

### Additional data
For some markdown documents, additional data may be required. See: `markdowns_for_analysis_and_figures/ASD_phenomics_figure_scripts/README.md`

### Configuration
If the data is all placed in the default location, then no changes to the configuration file should be required. Otherwise, make sure the `dataset` variable in `config.yml` points to where the raw data was extracted, or that the `data` variable points to where the processed data was downloaded (if the processed data was generated by running `main.R` then it doesn't need to be changed.)


<!-- USAGE EXAMPLES -->
## Usage
Raw data can be processed using `main.sh` (for `bash` environments) or `main.R`. Individual markdowns were generated using `knitr` in [RStudio](https://www.rstudio.com/). 

### Command line
The `bash` or `R` version of the main scripts do the same thing, except that `main.sh` also retains runtime logs in the `logs/` directory.

```
## Using bash
./main.sh

## Using Rscript
Rscript main.R
```
To run/re-run individual parts of the analysis, see the sub-scripts called in `main.R`.

<!-- LICENSE -->
<!--
## License

Distributed under the *** License. See `LICENSE` for more information.
-->


<!-- CONTACT -->
## Contact


Manuel Belmadani - [GitHub](https://github.com/mbelmadani) - [@Pragmatwit](https://twitter.com/pragmatwit) - manuel.belmadani@msl.ubc.ca  
Troy A. McDiarmid - [GitHub](https://github.com/troymcdiarmid) - [@Troy_McD_UBC](https://twitter.com/Troy_McD_UBC)

Further information and requests for resources and reagents should be directed to and will be
fulfilled by the Lead Contact, Catharine H. Rankin (crankin@psych.ubc.ca).


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
We would like to thank Dr. Evan L. Ardiel for useful comments and discussion regarding the manuscript. We would like to thank Dr. John Calarco, Dr. Erik Jorgensen, and Dr. Don Moerman and their labs for sharing their constructs and protocols or making them publicly available. We thank the C. elegans knockout consortium for several alleles, as well as the National Bioresource Project and the C. elegans Genetics Center (NIH Office of Research Infrastructure Programs, P40623 OD010440) for providing strains. We would like to thank Lexis Kepler and Anna Willms for assistance with behavioral and cloning experiments. We would like to thank Christine R. Ackerley for useful advice and discussions regarding figure design. We would also like to thank the Caenorhabditis Genetic Center and the National BioResource Project of Japan for strains. 

