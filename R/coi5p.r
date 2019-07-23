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






coi = function(str){
  #take a string as input and turin it into the decided upon data structure
    #store original str in $raw
    #build an aphid::DNAbin for the $raw sequence
}

set_frame = function(coi_obj){
  #input is the output structure from coi
  #set the reading frame and store the framed string in $framed
  #^this will require a lot of the code from the 5_correction_script to be called and for
  # the other functions in this package to be used (i.e. set_frame)

}

which_trans_table = function(str){
  #in translation.r - reads a string with the name of a taxonomic grouping and returns which
  #translation table should be used. 0 is returned to indicate censored translation when ambigious
}


translate = function(coi_obj, other_params){
  #take the string from the in_frame coi object and translate it to amino acids based on the
  #reading frame for the inframe objects and the translation table
  #store the output in the $AAseq for the COI object
}

check_coi = function(coi_obj){
  #takes the $AAseq and passes it through the aaPHMM, storing the output score.
  #this score can be used to give a probability that the sequence contains an error or not.
  #output save to $AAscore and $indel_likely (T/F) and $stop_codons (T/F)
}




