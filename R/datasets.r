###############################################################################
#' Nucleotide profile hidden Markov model for coi5p.
#'
#' This model is stored in the coi5p package and was trained on a representative
#' sample of the barcode of life database (http://www.boldsystems.org/index.php).
#'
#' @keywords internal
"nt_PHMM"
###############################################################################

###############################################################################
#' Amino acid profile hidden Markov model for coi5p.
#'
#' This model is stored in the coi5p package and was trained on a representative
#' sample of the barcode of life database (http://www.boldsystems.org/index.php).
#'
#' @keywords internal
"aa_PHMM"
###############################################################################

###############################################################################
#' Data frame containing the translation table recommendation.
#'
"trans_df"
###############################################################################

###############################################################################
#' Example coi5p DNA sequence string
#'
#' This string of barcode data is used in the package documentation's examples
#' and within the vignette demonstrating how to use the package.
#'
"example_nt_string"
###############################################################################

###############################################################################
#' Example barcode data.
#'
#' A nine line dataframe of coi5p barcode data with the following columns:
#'
#' id - the unique identifier for the sample
#'
#' genetic_code - the genetic code for translation of the sample (features NA for unknowns)
#'
#' taxa -  a taxonomic designation associated with the sample
#'
#' sequence - the DNA sequence associated with the sample
#'
#' notes - notes on the sequence structure
#'
"example_barcode_data"
###############################################################################
