# coi5p 
An R package for contextualization and evaluation of COI-5P barcode data
[![Build Status](https://travis-ci.com/CNuge/coi5p.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/coi5p)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
--------------------------------------------------------------------------

`coi5p` is an R package designed to aid users in the cleaning and analysis of COI-5P DNA barcode data.


## Installation

To download and install `coi5p` directly from cran, type the following commands into R.
```
install.packages("coi5p")
library(coi5p)
```
Alternatively, you can download the development version of `coi5p` directly from GitHub. You'll need to have the R package `devtools` installed and loaded.
```
#install.packages("devtools")
#library(devtools)
devtools::install_github("CNuge/coi5p", build_vignettes = TRUE)
library(coi5p)
```

## How to use it

Below is a brief demonstration to get the user started, please consult the package vignette for a more detailed explanation of `coi5p`'s functionality.

The package is built around the custom `coi5p` object, which takes a COI-5P DNA barcode sequence as input. The package contains functions for: 

  - setting the sequence in reading frame
  - translating the sequence to amino acids
  - checking the sequence for evidence of insertion or deletion errors

The basic `coi5p` analysis pipleine is as follows:
```
example_nt_string #an input DNA string, contained in the coi5p package for demonstration purposes

#step 1: build the coi5p object
dat = coi5p(example_nt_string, name="example_sequence_1")

#step 2: frame the sequence
dat = frame(dat)

#step 3: by default censored translation is performed - see vignette for details
dat = translate(dat)

##step 3a: if taxonomy is known, but the translation table is not, a helper function
#can be used to look up the proper translation table.
which_trans_table("Scyliorhinidae")

#step 3a: the proper transaltion table can be passed to the translation function
dat = translate(dat, trans_table = 2)

#step 4: check to see if an insertion or deletion is likely
dat = indel_check(dat)
dat
```
All of the steps of the pipeline can be called at once through the `coi5p_pipe` function.
```
output = coi5p_pipe(example_nt_string)
output
```
Afterwards, the coi5p object has the following components that can be called by the user.
```
output$name         #the name of the sequence 
output$raw          #the input DNA sequence
output$framed       #the DNA sequence set in reading frame
output$aaSeq        #the amino acid sequence
output$aaScore      #the log likelihood score of the amino acid sequence - see vignette for details
output$indel_likely #a boolean indicating whether the sequence should be double checked for indel errors
output$stop_codons  #a boolean indicating whether the amino acid sequence contains stop codons.
output$data         #used internally by the function for data storage
```
Most use cases will involve the analysis of multiple sequences. Please consult the package's vignette for a suggested workflow for batch analysis and demonstration of how the batch analysis helper function can be used to build dataframes out of multiple coi5p objects.
