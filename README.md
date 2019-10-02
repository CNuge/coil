# coil 
An R package for contextualization and evaluation of COI-5P barcode data
[![Build Status](https://travis-ci.com/CNuge/coi5p.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/coi5p)
[![codecov](https://codecov.io/gh/CNuge/coi5p/branch/master/graph/badge.svg)](https://codecov.io/gh/CNuge/coi5p)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
--------------------------------------------------------------------------

**coil** is an R package designed for the cleaning, contextualization and assessment of cytochrome c oxidase I DNA barcode data ([COI-5P, or the five prime portion of COI](https://en.wikipedia.org/wiki/Cytochrome_c_oxidase_subunit_I)). It contains functions for placing COI-5P barcode sequences into a common reading frame, translating DNA sequences to amino acids and for assessing the likelihood that a given barcode sequence includes an insertion or deletion error. These functions are provided as a single function analysis pipeline and are also available individually for efficient and targeted analysis of barcode data.

## Installation

At the moment, you can download the development version of `coil` directly from GitHub. You'll need to have the R package `devtools` installed and loaded. Also note if the `build_vignettes` option is set to true, you will need to have the R package `knitr` installed.
```
#install.packages("devtools")
#install.packages("knitr") #required if build_vignettes = TRUE
#library(devtools) 
devtools::install_github("CNuge/coil", build_vignettes = TRUE)
library(coil)
```

The vignette can then be accessed from R using the following command:
```
vignette("coil-vignette")
```

## How to use it

Below is a brief demonstration to get the user started, please consult the package vignette for a more detailed explanation of `coil`'s functionality.

The package is built around the custom `coi5p` object, which takes a COI-5P DNA barcode sequence as input. The package contains functions for: 

  - setting the sequence in reading frame
  - translating the sequence to amino acids
  - checking the sequence for evidence of insertion or deletion errors

The basic `coi5p` analysis pipeline is as follows:
```
example_nt_string #an input DNA string, contained in the coil package for demonstration purposes

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
```
Calling the variable name prints the coi5p object's summary and shows all of the important information, including: the original raw sequence, the sequence set in reading frame, the amino acid sequence and the summary stats regarding the likelihood of the sequence containing an error.
```
output 
#calling output will return the following:
#coi5p barcode sequence
#raw sequence:
#ctctacttgatttttggtgcatgag...ggacccaattctctatcaacactta
#framed sequence:
#---ctctacttgatttttggtgcat...ggacccaattctctatcaacactta
#Amino acid sequence:
#-LYLIFGAWAG?VG?ALSLLIRAEL...LTDRNLNTTFFDPAGGGDPILYQHL
#The sequence likely does not contain an insertion or deletion.
#Stop codon present: FALSE, Amino acid PHMM score:-206.22045
```
The coi5p object has the following components that can be extracted by the user using the dollar sign notation.
```
output$name         #the name of the sequence 
output$raw          #the input DNA sequence
output$framed       #the DNA sequence set in reading frame
output$aaSeq        #the amino acid sequence
output$aaScore      #the log likelihood score of the amino acid sequence - see vignette for details
output$indel_likely #a boolean indicating whether the sequence should be double checked for indel errors
output$stop_codons  #a boolean indicating whether the amino acid sequence contains stop codons.
output$data         #contains the generated nucleotide and amino acid hidden state paths.
```
Most use cases will involve the analysis of multiple sequences. Please consult the package's vignette for a suggested workflow for batch analysis and demonstration of how the batch analysis helper function can be used to build dataframes out of multiple coi5p objects.
