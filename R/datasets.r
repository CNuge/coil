
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

#was accessing like so:
# coi5p:::nt_PHMM
# but this throws a note on compilation of the package

###############################################################################
#' Nucleotide profile hidden markov model for coi5p.
#'
#' removed param: phmm the profile hidden markov model against which the coi5p class should
#' be framed. By defualt the function will use the nt_PHMM variable
#' stored in the coi5p package, which was trained on a representitive sample of the
#' barcode of life database (). A user may wish to use a custom derived PHMM, in which
#' case they should consult the aphid package () for custom PHMM derivation.
#'
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
"aa_PHMM"
###############################################################################

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



###
# For storing data in the data folder -fallback
###

# save(trans_df, file = 'data/trans_df.RData')
# save(nt_PHMM, file = 'data/nt_PHMM.RData')
# save(aa_PHMM, file = 'data/aa_PHMM.RData')
#
#globalVariables(c("nt_PHMM", "aa_PHMM", "translation_table_data"))
