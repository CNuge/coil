
###############################################################################
#' Nucleotide profile hidden markov model for coi5p.
#'
#' removed param: phmm the profile hidden markov model against which the coi5p class should
#' be framed. By defualt the function will use the nt_PHMM variable
#' stored in the coi5p package, which was trained on a representitive sample of the
#' barcode of life database (). A user may wish to use a custom derived PHMM, in which
#' case they should consult the aphid package () for custom PHMM derivation.
#'
#' @keywords internal
"nt_PHMM"
###############################################################################

###############################################################################
#' Amino acid profile hidden markov model for coi5p.
#'
#' removed param: phmm the amino acid profile hidden markov model against which the coi5p class
#' should be checked for errors. By defualt the function will use the aa_PHMM variable
#' stored in the coi5p package, which was trained on a representitive sample of the
#' barcode of life database (). A user may wish to use a custom derived PHMM, in which
#' case they should consult the aphid package () for custom PHMM derivation.
#' @keywords internal
"aa_PHMM"
###############################################################################


#TODO - found one error in the data from Tyler - Pyuridae had wrong translation table
#should be 13 but was 11... need to do a big double check and make sure there aren't more erroneous ones.
#then change here if needed and resave to the file
#example below done for Pyuridae
# load('R/sysdata.rda')
# trans_df$trans_table[trans_df$taxon == 'Pyuridae'] = 13
# trans_df$trans_table[trans_df$taxon ==	'Stolidobranchia'] = 0

# added this to it:
# trans_df$trans_table[trans_df$taxon ==	'Ascidiacea'] = 0
#asc = data.frame(trans_table = 0, taxon = "Ascidiacea", level = "class")
#trans_df = rbind(trans_df, asc)
# head(trans_df)
#need to look for more missing information!
###############################################################################
#' Data frame containing the translation table recommendation.
#'
"trans_df"
###############################################################################

###############################################################################
#' Example coi5p DNA sequence string
#' example_nt_string = 'ctctacttgatttttggtgcatgagcaggaatagttggaatagctttaagtttactaattcgcgctgaactaggtcaacccggatctcttttaggggatgatcagatttataatgtgatcgtaaccgcccatgcctttgtaataatcttttttatggttatacctgtaataattggtggctttggcaattgacttgttcctttaataattggtgcaccagatatagcattccctcgaataaataatataagtttctggcttcttcctccttcgttcttacttctcctggcctccgcaggagtagaagctggagcaggaaccggatgaactgtatatcctcctttagcaggtaatttagcacatgctggcccctctgttgatttagccatcttttcccttcatttggccggtatctcatcaattttagcctctattaattttattacaactattattaatataaaacccccaactatttctcaatatcaaacaccattatttgtttgatctattcttatcaccactgttcttctactccttgctctccctgttcttgcagccggaattacaatattattaacagaccgcaacctcaacactacattctttgaccccgcagggggaggggacccaattctctatcaacactta'
"example_nt_string"
###############################################################################


###############################################################################
#' Example barcode data.
#' A nine line dataframe with the following columns:
#' id - the unique identifier for the sample
#' genetic_code - the genetic code for translation of the sample (features NA for unknowns)
#' taxa -  a taxonomic designation associated with the sample
#' sequence - the DNA sequence associated with the sample
#' notes - notes on the sequence structure
#'
"example_barcode_data"
###############################################################################


# For dev/testing purposes only
# library(ape)
# library(aphid)
# library(seqinr)
# source('R/deploy_PHMMs.r')
# source('R/translation.r')
# load('R/sysdata.rda')

###
# For storing data in the R folder - used currently
###
#
#use_data(nt_PHMM , aa_PHMM, trans_df, example_nt_string, overwrite = TRUE, internal = TRUE)

#This saves the data to the /data/ folder so it can be accessed by the package
#seems better to put them in the sysdata.rda folder and use them internall

