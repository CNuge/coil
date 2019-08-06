#' Run the entire coi5p pipeline for am input sequence.
#'
#' This function will take a raw DNA sequence string and run each of the coi5p methods in turn
#' (coi5p, frame, translate, indel_check). Note that if you are not intersted in all components
#' of the output (i.e. only want sequences in frame or translated), then the coi5p analysis functions
#' can be called individually to avoid unnecessary computation.
#'
#'
#' @param x a nucleotide string.
#' Valid characters within the nucleotide string are: a,t,g,c,-,n.
#' The nucleotide string can be input as upper case, but will be automatically converted to lower case.
#' @param ... additional arguments to be passed between methods.
#' @param name an optional character string. Identifier for the sequence.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids.
#' Default is 0, which indicates that censored translation should be performed. If the taxonomy
#' of the sample is known, use the function which_trans_table() to determine the translation table to use.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would make the second character in the framed sequence the the first nucelotide of
#' the first codon.
#' @param indel_threshold the log likelihood threshold used to assess whether or not sequences
#' are likely to contain an indel. Default is -345.95. Values lower than this will be classified
#' as likely to contain an indel and values higer will be classified as not likely to contain an indel.
#'
#' @return an object of class code{"coi5p"}
#' @examples
#' dat = coi5p_pipe(example_nt_string )
#' #full coi5p object can then be printed
#' dat
#' #components of output coi5p object can be called individually:
#' dat$raw    #raw input sequence
#' dat$name   #name that was passed
#' dat$framed #sequence in common reading frame
#' dat$aaSeq  #sequence translated to amino acids (censored)
#' dat$indel_likely #whether an insertion or deletion likely exists in the sequence
#' dat$stop_codons #whether or not there are stop codons in the amino acid sequence
#'
#' dat = coi5p_pipe(example_nt_string , trans_table = 5)
#' dat$aaSeq #sequence translated to amino acids using designated translation table
#' @seealso \code{\link{coi5p}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{translate}}
#' @seealso \code{\link{indel_check}}
#' @seealso \code{\link{which_trans_table}}
#' @name coi5p_pipe
#' @export
coi5p_pipe = function(x, ... ,
                      name = character(),
                      trans_table = 0,
                      frame_offset = 0,
                      indel_threshold = -346.95){
  dat = coi5p(x, name=name)
  dat = frame(dat )
  dat = translate(dat, trans_table = trans_table, frame_offset = frame_offset)
  dat = indel_check(dat, indel_threshold=indel_threshold)
  return(dat)
}
