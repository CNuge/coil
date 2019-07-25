# move the code from the corrections script in here and structure it to work on
# a standardized object

#possibly this for getting the data for the PHMMs into the proper fmt
#.onLoad() = function(libname, pkgname)

#start with r/ chapter of textbook and begin work accordinlgy. Go through R studio

#use: package::function() syntax for external functions so as to make it explicicit that that is needed

#NOTE FOR GETTING PHMMS and TRANS TABLES in:
# aphid has a folder /data with .RData (binary?) objects within that store the data, presumably I'll need ot do this.


#' This is how to document in a way that can be accessed by ?coi
#' Add roxygen comments to your .R files.
#' Run devtools::document() (or press Ctrl/Cmd + Shift + D in RStudio) to convert roxygen comments to .Rd files. (devtools::document() calls roxygen2::roxygenise() to do the hard work.)
#' Preview documentation with ?.
#' Rinse and repeat until the documentation looks the way you want.
#' It will generate a man/foo.Rd file, which should not be modified by hand
#' There is shorthand for arguments, examples etc. see the book, most common for functions:
#' @param
#' @examples
#' @return
foo = function(x){

}


# building the data and functions I've created for manipulating COI-5P sequences
# into a generic s3 function

# For dev/testing purposes only
# library(ape)
# library(aphid)
# library(seqinr)
# source('R/deploy_PHMMs.r')
# source('R/translation.r')
# nt_phmm_file = './required_data/COI5P_nt.PHMM'
# aa_phmm_file = './required_data/COI5P_aa.PHMM'
# trans_df = read.table('./required_data/family_tanslation_table.tsv', header = TRUE, sep = '\t')
# nt_PHMM = readPHMM(nt_phmm_file)
# aa_PHMM = readPHMM(aa_phmm_file)

#This saves the data to the /data/ folder so it can be accessed by the package
# save(trans_df, file = 'data/trans_df.RData')
# save(nt_PHMM, file = 'data/nt_PHMM.RData')
# save(aa_PHMM, file = 'data/aa_PHMM.RData')
#
# Generating the namespace with roxygen2 is just like generating function documentation with roxygen2.
# You use roxygen2 blocks (starting with #') and tags (starting with @).
#  The workflow is the same:
#   Add roxygen comments to your .R files.
#   Run devtools::document() (or press Ctrl/Cmd + Shift + D in RStudio) to convert roxygen comments to .Rd files.
#   Look at NAMESPACE and run tests to check that the specification is correct.
#   Rinse and repeat until the correct functions are exported.

#To export an object, put @export in its roxygen block - just don't do this for the functions the user doesn't have to see.
# ^not needed for data, these should just be avaliable?
# TODO - document functions and make sure only the user facing ones are exported
# TODO - take the positions where functions from other libraries are used, use them in the tidyverse::func() style


########################
# coi5p - Initialization of the class


# three functions should be provided at minimum"

# constructor - efficiently creates new objects with the correct structure - new_coi5p
# validator - perform computationally expensive checks to make sure the obj has the correct vals - validate_coi5p
# helper - provide a way for other to create objects of the class - coi5p

new_coi5p = function(x = character(), name = character()){
  stopifnot(is.character(x))

  structure(list(name = name, raw = tolower(x)) , class = "coi5p")
}


# take a new instance and run validation checks on the sequence
# make sure the sequence has only ATGCN-
# make sure the sequence has length greater than zero
validate_coi5p = function(new_instance){

  #TODO - check that the string is composed of the correct character types
  new_instance
}


#'
#'@param x a nucleotide string.
#'@param  name an optional character string. Identifier for the sequence.
#'
#'@return an object of class code{"coi5p"}
#'@examples
#'
#'@export
coi5p = function(x = character(), name = character()){

  # TODO - coerce the input into a lower case character string
  # TODO - vector of characters can be another acceptable input

  #if vector, paste them together
  validate_coi5p(new_coi5p(x, name))

}




###########################
# coi5p - Generics and methods

# TODO - add checks to make sure data structures required have been initialized
# if not then return a warning saying that the previous method needs to be run first


#' ! this is where we would document this function in detail
#'@param x
#'
#'@return an object of class code{"coi5p"}
#'@examples
#'
#'
#'@export
frame = function(x){
  UseMethod("frame")
}


frame.coi5p = function(x){
  #input is the output structure from coi
  #set the reading frame and store the framed string in $framed

  x$data$ntBin = individual_DNAbin(x$raw)
  x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)

  if(leading_ins(x$data$ntPHMMout[['path']])){
    trim_temp  = set_frame(x$raw, x$data$ntPHMMout[['path']])
    x$data$ntBin = individual_DNAbin(trim_temp)
    x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
  }else{
    trim_temp = x$raw
  }

  x$framed = set_frame(trim_temp, x$data$ntPHMMout[['path']])

  return(x)
}

#'
#'
#'@param x a coi5p class object for which frame() has been run
#'@param trans_table
#'
#'@return an object of class code{"coi5p"}
#'@examples
#'
#'@export
translate = function(x, ...){
  UseMethod("translate")
}


translate.coi5p = function(x, trans_table = 0){
  if(trans_table == 0){
    x$aaSeq = censored_translation(x$framed)
  }else{
    #split the DNA string into a vector, all characters to lower case
    dna_list = strsplit(gsub('-', 'n', as.character(tolower(x$framed))),"")
    dna_vec = dna_list[[1]]
    #translate using the designated numcode, returns a vector of AAs
    aa_vec = seqinr::translate(dna_vec, frame = 0, numcode=trans_table, ambiguous= TRUE, NAstring = '-')

    x$aaSeq = paste(aa_vec, collapse= "")
  }

  return(x)
}


#' Check a translated coi5p sequence to see if an indel error is likely present
#'
#'
#'@param x a coi5p class object for which frame() and translate() have been run
#'@param indel_threshold
#'
#'@return an object of class code{"coi5p"}
#'@examples
#'
#'@export
indel_check = function(x, ...){
  UseMethod("indel_check")
}

indel_check.coi5p = function(x, indel_threshold = -346.95 ){

  x$data$aaBin = individual_AAbin(x$aaSeq)
  x$data$aaPHMMout = aphid::Viterbi(aa_PHMM, x$data$aaBin, odds = FALSE)

  x$AAscore = x$data$aaPHMMout[['score']] #have this print a threshold
  if(x$AAscore > indel_threshold){
    x$indel_likely = FALSE
  }else{
    x$indel_likely = TRUE
  }

  if(grepl('\\*', x$aaSeq)){
    x$stop_codons = TRUE
  }else{
    x$stop_codons = FALSE
  }

  return(x)
}

