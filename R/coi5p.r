# move the code from the corrections script in here and structure it to work on
# a standardized object

#possibly this for getting the data for the PHMMs into the proper fmt
#.onLoad() = function(libname, pkgname)

#start with r/ chapter of textbook and begin work accordinlgy. Go through R studio

#use: package::function() syntax for external functions so as to make it explicicit that that is needed

#NOTE FOR GETTING PHMMS and TRANS TABLES in:
# aphid has a folder /data with .RData (binary?) objects within that store the data, presumably I'll need ot do this.

#TODO - run examples with N in sequence, confirm that these are being left alone as unknonwns/placeholders
# that are treated differently than the -
# if this is not the case then modify the wrapper functions for aphid accordingly



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

#data structure - use an s3 or s4 class to store the data associated with the functions in a standardized way
#may be a good idea to associate all of the functions with the class, have the PHMMs be stored within the class?
#

#^ read these two sections of adv-r and decide which is more appropriate
# s3 seems very flexible and recommended
# s4 has documented slots so may be better to use that to store the data in its different forms




# building the data and functions I've created for manipulating COI-5P sequences
# into a generic s3 function

# For testing purposes only
# library(ape)
# library(aphid)
# library(seqinr)
#   source('deploy_PHMMs.r')
# source('translation.r')
#   nt_phmm_file = '../required_data/COI5P_nt.PHMM'
# aa_phmm_file = '../required_data/COI5P_aa.PHMM'
# nt_PHMM = readPHMM(nt_phmm_file)
# aa_PHMM = readPHMM(aa_phmm_file)



########################
# coi5p - Initialization of the class


# three functions should be provided at minimum"

# constructor - efficiently creates new objects with the correct structure

# validator - perform computationally expensive checks to make sure the obj has the correct vals

#helper - provide a way for other to create objects of the class


new_coi5p = function(str = character()){
  stopifnot(is.character(str))

  structure(list(raw = tolower(str)) , class = "coi5p")
}


# take a new instance and run validation checks on the sequence
# make sure the sequence has only ATGCN-
# make sure the sequence has length greater than zero
validate_coi5p = function(new_instance){


  new_instance
}



#helper function - this is the user facing part of the class
coi5p = function(str = character()){
  
  #coerce the input into a lower case character string
  # vector can be another acceptable input

  #if vector, paste them together 
  validate_coi5p(new_coi5p(str))

}




###########################
# coi5p - Generics and methods

# TODO - add checks to make sure data structures required have been initialized
# if not then return a warning saying that the previous method needs to be run first


#' this is where we would document this function in detail
#' this is a method dispatch, its job is to find the specific 
#' implementation for this class
frame = function(coi_obj){
  UseMethod("frame")
}


#' this method version calls ins_front_trim from the deploy PHMMs file, can that function be
#' hidden from the user so it only exists in the background?
frame.coi5p = function(coi_obj){
  #input is the output structure from coi
  #set the reading frame and store the framed string in $framed
  #^this will require a lot of the code from the 5_correction_script to be called and for
  # the other functions in this package to be used (i.e. set_frame)


#######
# should these be separate function stored elsewhere? 
# Or are nested functions within the method fine?
# 
  # TODO - run the individual_DNAbin() function to build the class's $DNAbin
  coi_obj$ntBin = individual_DNAbin(coi_obj$raw)
  # TODO - pass the $DNAbin into Viterbi and store as ntPHMMout
  # TODO: !!!!!!!FIGURE OUT WHERE nt_PHMM WILL BE STORED!!!!!!!!!
  coi_obj$ntPHMMout = Viterbi(nt_PHMM, coi_obj$ntBin, odds = FALSE)

  if(leading_ins(coi_obj$ntPHMMout[['path']])){
    trim_temp  = set_frame(coi_obj$raw, coi_obj$ntPHMMout[['path']])
    coi_obj$ntBin = individual_DNAbin(trim_temp)
    coi_obj$ntPHMMout = Viterbi(nt_PHMM, coi_obj$ntBin, odds = FALSE)
  }else{
    trim_temp = coi_obj$raw
  }

  coi_obj$framed = set_frame(trim_temp, coi_obj$ntPHMMout[['path']])

  return(coi_obj)
}

#`
#`
#`
#`
#`
translate = function(coi_obj, ...){
  UseMethod("translate")
}


translate.coi5p = function(coi_obj, trans_table = 0){
  if(trans_table == 0){
    coi_obj$aaSeq = censored_translation(coi_obj$framed)
  }else{
    #split the DNA string into a vector, all characters to lower case
    dna_list = strsplit(gsub('-', 'n', as.character(tolower(coi_obj$framed))),"")
    dna_vec = dna_list[[1]]
    #translate using the designated numcode, returns a vector of AAs
    aa_vec = seqinr::translate(dna_vec, frame = 0, numcode=trans_table, ambiguous= TRUE, NAstring = '-')

    coi_obj$aaSeq = paste(aa_vec, collapse= "")
  }

  return(coi_obj)
}


#' have a default indel threshold learned from the real world data
#' checks for indels likelihood, reports the path score and also checks for stop codons
indel_check = function(coi_obj){
  UseMethod("indel_check")
}

indel_check.coi5p = function(coi_obj, indel_threshold = -346.95 ){
  
  coi_obj$aaBin = individual_AAbin(coi_obj$aaSeq)
  coi_obj$aaPHMMout = Viterbi(aa_PHMM, coi_obj$aaBin, odds = FALSE)

  coi_obj$AAscore = coi_obj$aaPHMMout[['score']] #have this print a threshold
  if(coi_obj$AAscore > indel_threshold){
    coi_obj$indel_likely = FALSE
  }else{
    coi_obj$indel_likely = TRUE   
  }

  if(grepl('\\*', coi_obj$aaSeq)){
    coi_obj$stop_codons = TRUE
  }else{
    coi_obj$stop_codons = FALSE
  }

  return(coi_obj)
}


#TODO:
# Write some unit tests for it here using input COI sequences:
  # a short one
  #one with an in
  #one with a del
  #one with lots of leading bp

# Once the unit tests are done and I can prove this class setup is working, then
# move it over to coi5p and begin documenting and cleaning up the code for 
# package development


seq_normal = 'cttcacttgatttttggtgcatgagcaggaatagtaggaactgctttaagtctccttattcgagcagaactgggtcaacctggttcacttttaggtgatgaccagatctacaatgtgatcgtaaccgcccatgctttagtaataattttttttatagttataccggtaataattggtggctttggaaactgactagtgcccctaataattggtgcaccagatatggcctttcctcgaataaataacataagtttttgactccttccaccatcattccttttattattagcttctgcaggggtagaagccggagctggcaccggctgaacagtttacccacccttatcgggtaatttagcacatgccgggccatctgttgatttaactattttttcacttcatttagcaggtgtatcatcaattttagcctcaattaattttatcacaactattattaatataaaaccaccagctatttctcaataccaaacaccattatttgtttgatccgttcttgtaactactattttactacttttagcccttccagtacttgcagctggaattacaatattattaacagatcgaaacctaaataccacattctttgaccctgctggtggaggagatcctatccactatcaacatcta'
seq_short = 'CCAGGTCTATAACGTAGTCGTCACAGCCCATGCCTTCGTAATAATCTTCTTCATAGTTATGCCTATTATAATCGGAGGATTTGGGAACTGACTAGTCCCTCTAATAATCGGAGCCCCAGACATAGCATTCCCACGAATAAACAACATAAGCTTCTGACTACTCCCCCCATCGTTCCTCCTACTACTAGCGTCCTCTACTGTAGAAGCAGGAGTTGGCACAGGATGAACAGTATACCCACCATTAGCCGGCAACTTAGCCCACGCTGGAGCTTCAGTTGACTTAGCAATCTTCT'
seq_in = 'ctttatttaatttttggtgcatgagcaggaatagttggaacggctttaagtcttctaatccgagctgaactaggacaacctgggtctctcctaggggatgatcaaatttataatgtaattgtaaccgcccatgcttttgtaataattttctttatagtaatacctgtcataattggtggttttggaaattaactaattccattaataattggtgcacctgacatagccttcccacgaataaataacataagctcctgacttcttccaccatcatttctccttctcctcgcctccgctggggttgaagccggagcaggtaccggttgaacagtttaccccccactggcaagcaaccttgctcatgccggaccatctgttgatttagctatcttctccctccatttagctggtatttcatcaattttagcctcaatccaacttcatcacaactattattaatataaaacccccagccatttctcaatatcaaacaccactatttgtttgatctatccttgtaactactattcttctcctcctttccctcccagttcttgcagcaggaattacaatcttacttacagaccgcaaccttaatactacattctttgatcctgcaggtggaggagacccaatcctttaccaacaccta'
seq_del = 'ctttacttaatctttggtgcatgagcaggaatagtaggtacagcccttagcttgcttattcgagcagaattaagccaacctggcacactcctgggagacgatcagatctacaatgttatcgtaactgctcacgcttttgtaataattttttttatggttatacctgtaataattggtgggttcggaaactgattagtgcctttaataattggtgcaccggacatagctttcccacgaataaataacataagcttttgactgctacccccctccctcctattacttttggcctctgctggagttgaagccggagccggaactggttgaacagtttatccccccctcgcaagtaatatagcccacgctggggcatcagtagacttagctattttctcgctccatttagcggtatttcctcaattcttgcctctatcaactttattacaaccattattaatataaaaccgcctgccatctctcaatatcaaacacccctttttgtttgatctattcttgtaaccacagtcctactcctcctttcacttcctgttcttgcagccgcaattacaatactacttaccgaccgtaatttaaacacaacattttttgatcctgctggtgggggtgacccaattctttaccaacattta'
seq_long = 'AACCGCTGATTATTTTCAACCAACCACAAAGATATCGGCACCCTTTACCTTCTATTTGGTGCCTGAGCTGGTATAGTAGGAACCGCCTTAAGCCTACTAATTCGCGCCGAACTAGGCCAACCCGGAACTCTACTCGGAGATGACCAAATCTACAACGTAATTGTAACCGCACATGCATTTGTAATAATTTTCTTTATAGTAATGCCTATTATAATCGGTGGATTCGGCAACTGACTAGTTCCTCTGATAATTGGAGCCCCTGATATAGCATTTCCTCGGATAAATAACATAAGCTTTTGACTTCTTCCCCCATCTTTCCTGTTACTCCTAGCATCCTCTATGGTTGAGGCCGGAGCAGGAACAGGTTGTAAGGTAGTATTCTCTCTAGCAGGCAACCTAGCCCATGCAGGAGCCTCAGTAGACCTAACTATTTTCTCCCTACACCTGGCAGGTGTCTCTTCAATTCTAGGAGCCATTAATTTTATTACAACTATTATTAATATAAAACCCCCTGCGATGTCACAGTATCAAACCCCCTTGTTTGTATGATCTGTACTAATCACTGCCGTACTTCTCCTTCTCTCACTTCCTGTATTAGCAGCTGGTATCACAATACTACTAACGGACCGAAACCTGAACACAACCTTTTTTGACCCAGCAGGAGGAGGAGACCCTATCCTATATCAACACCTATTCTGATTCTTTGGGCACCCTGAAGTATATATTCTTATTTTACCTGGGTTTGGGATAATCTCCCATATTGTGACCTACTATTCAGGAAAAAAAGAACCATTCGGATATATAGGAATAGTATGAGCCATAATATCAATTGGGTTCCTAGGATTCATTGTATGAGCCCACCATATATTCACAGTCGGAATAG'


test_normal = coi5p(seq_normal)
test_normal
test_normal = frame(test_normal)
test_normal
test_normal = translate(test_normal)
test_normal
test_normal = translate(test_normal, trans_table = 5)
test_normal
test_normal = indel_check(test_normal)
test_normal


test_short = coi5p(seq_short)
test_short
test_short = frame(test_short)
test_short
test_short = translate(test_short)
test_short
test_short = translate(test_short, trans_table = 5)
test_short
test_short = indel_check(test_short)
test_short


test_in = coi5p(seq_in)
test_in
test_in = frame(test_in)
test_in
test_in = translate(test_in)
test_in
test_in = translate(test_in, trans_table = 5)
test_in
test_in = indel_check(test_in)
test_in



test_del = coi5p(seq_del)
test_del
test_del = frame(test_del)
test_del
test_del = translate(test_del)
test_del
test_del = translate(test_del, trans_table = 5)
test_del
test_del = indel_check(test_del)
test_del



test_long = coi5p(seq_long)
test_long
test_long = frame(test_long)
test_long
test_long = translate(test_long)
test_long
test_long = translate(test_long, trans_table = 5)
test_long
test_long = inlong_check(test_long)
test_long

