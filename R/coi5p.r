########################
# coi5p - Initialization of the class

#' Build a new coi5p class instance.
#'
#' @keywords internal
new_coi5p = function(x = character(), name = character()){
  stopifnot(is.character(x))
  stopifnot(is.character(name))
  if(length(x) == 0){
    stop("Must pass a DNA sequence.")
  }

  structure(list(name = name, raw = tolower(x)) , class = "coi5p")
}

#' Validate the new coi5p class instance.
#'
#' @keywords internal
validate_coi5p = function(new_instance){
  # take a new instance and run validation checks on the sequence
  # make sure the sequence has only ATGCN-
  # make sure the sequence has length greater than zero
  allowed = c("-", "a", "c", "g", "n","t")
  for(c in sort(unique(strsplit(new_instance$raw, "")[[1]]))){
    if(!c %in% allowed){
      stop(paste("Unallowed character in DNA string:", c,
                 "\n Valid characters are: a t g c - n"))
    }
  }

  if(new_instance$raw != tolower(new_instance$raw)){
    stop("Unable to convert DNA string to lower case")
  }

  new_instance
}


#' Build a coi5p object from a DNA sequence string.
#'
#' @param x a nucleotide string.
#' Valid characters within the nucleotide string are: a,t,g,c,-,n.
#' The nucleotide string can be input as upper case, but will be automatically converted to lower case.
#' @param name an optional character string. Identifier for the sequence.
#'
#' @return an object of class code{"coi5p"}
#' @examples
#' dat = coi5p(example_nt_string)
#' #named coi5p sequence
#' dat = coi5p(example_nt_string, name = "example_seq1")
#' #components in output coi5p object:
#' dat$raw
#' dat$name
#' @name coi5p
#' @export
coi5p = function(x = character(), name = character()){
  validate_coi5p(new_coi5p(tolower(x), name))
}



###########################
# coi5p - Generics and methods

#' Take a coi5p sequence and place it in reading frame.
#'
#' @param x a coi5p class object
#' @param ... additional arguments to be passed between methods.
#'
#' @return an object of class code{"coi5p"}
#' @seealso \code{\link{coi5p}}
#' @examples
#' #previously run function:
#' dat = coi5p(example_nt_string )
#'
#' dat = frame(dat)
#'
#' #additional components in output coi5p object:
#' dat$framed
#' @export
#' @name frame
frame = function(x, ...){
  UseMethod("frame")
}

####
#' @rdname frame
#' @export
frame.coi5p = function(x, ... ){
  #input is the output structure from coi5p
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

#' Translate a coi5p sequence.
#'
#' @param x a coi5p class object for which frame() has been run.
#' @param ... additional arguments to be passed between methods.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids.
#' Default is 0, which indicates that censored translation should be performed. If the taxonomy
#' of the sample is known, use the function which_trans_table() to determine the translation table to use.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would make the second character in the framed sequence the the first nucelotide of
#' the first codon.
#'
#' @return an object of class code{"coi5p"}
#' @seealso \code{\link{coi5p}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{which_trans_table}}
#' @examples
#' #previously run functions:
#' dat = coi5p(example_nt_string )
#' dat = frame(dat)
#' #translate when the translation table is not known:
#' dat = translate(dat)
#' #translate when the translation table is known:
#' dat = translate(dat, trans_table = 5)
#' #additional components in output coi5p object:
#' dat$aaSeq
#'@name translate
translate = function(x, ...){
  UseMethod("translate")
}

####
#' @rdname translate
#' @export
translate.coi5p = function(x, ..., trans_table = 0, frame_offset = 0){
  if(is.null(x$framed)){
    stop("translate function only accepts framed coi5p objects. See function: frame.")
  }

  if(trans_table == 0){
    x$aaSeq = censored_translation(x$framed, reading_frame = (frame_offset+1))
  }else{
    #split the DNA string into a vector, all characters to lower case
    dna_list = strsplit(gsub('-', 'n', as.character(tolower(x$framed))),"")
    dna_vec = dna_list[[1]]
    #translate using the designated numcode, returns a vector of AAs
    aa_vec = seqinr::translate(dna_vec, frame = frame_offset, numcode=trans_table, ambiguous= TRUE, NAstring = '-')

    x$aaSeq = paste(aa_vec, collapse= "")
  }

  return(x)
}


#' Check if coi5p sequence likely contains an indel error.
#'
#' @param x a coi5p class object for which frame() and translate() have been run.
#' @param ... additional arguments to be passed between methods.
#' @param indel_threshold the log likelihood threshold used to assess whether or not sequences
#' are likely to contain an indel. Default is -363.87. Values lower than this will be classified
#' as likely to contain an indel and values higer will be classified as not likely to contain an indel.
#'
#' @return an object of class code{"coi5p"}
#' @seealso \code{\link{coi5p}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{translate}}
#' @examples
#' #previously run functions:
#' dat = coi5p(example_nt_string)
#' dat = frame(dat)
#' dat = translate(dat)
#' #current function
#' dat = indel_check(dat)
#' #with custom indel threshold
#' dat = indel_check(dat, indel_threshold = -400)
#' #additional components in output coi5p object:
#' dat$stop_codons #Boolean - Indicates if there are stop codons in the amino acid sequence.
#' dat$indel_likely #Boolean - Indicates if there is likely a insertion or deletion in the sequence.
#' dat$aaScore #view the amino acid log likelihood score
#' @name indel_check
indel_check = function(x, ...){
  UseMethod("indel_check")
}

####
#' @rdname indel_check
#' @export
indel_check.coi5p = function(x, ..., indel_threshold = -363.87){
  if(is.null(x$framed)|is.null(x$aaSeq) ){
    stop("indel_check function accepts framed and translated coi5p objects. See functions: frame, translate.")
  }

  x$data$aaBin = individual_AAbin(x$aaSeq)
  x$data$aaPHMMout = aphid::Viterbi(aa_PHMM, x$data$aaBin, odds = FALSE)

  x$aaScore = x$data$aaPHMMout[['score']]
  if(x$aaScore > indel_threshold){
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
