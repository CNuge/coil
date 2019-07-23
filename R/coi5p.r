# move the code from the corrections script in here and structure it to work on
# a standardized object

#possibly this for getting the data for the PHMMs into the proper fmt
#.onLoad() = function(libname, pkgname)

#start with r/ chapter of textbook and begin work accordinlgy. Go through R studio

#use: package::function() syntax for external functions so as to make it explicicit that that is needed

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




