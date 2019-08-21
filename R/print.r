#' print summary of coi5p object.
#' @param x a coi5p class object.
#' @param ... additional arguments to be passed between methods.
#' @keywords internal
print.coi5p = function(x, ...){

  l1 = "coi5p barcode sequence"

  #addition if name
  if(length(x$name) != 0 ){
    l2 = paste(": ", x$name, "\n",
               sep = "")
  }else{
    l2 = "\n"
  }

  l3 = paste("raw sequence:\n",
              substr(x$raw,1, 25), "..." , substr(x$raw, (nchar(x$raw)-24), nchar(x$raw)),
             "\n",
              sep = "")

  lines = c(l1, l2, l3)

  #addition if framed
  if( "framed" %in% names(x) ){
    l4 = paste("framed sequence:\n",
               substr(x$framed,1, 25), "...", substr(x$framed, (nchar(x$framed)-24), nchar(x$framed)),
               "\n",
               sep = "")
    lines = c(lines, l4)
  }

  #addition if translated
  if( "aaSeq" %in% names(x) ){
    l5 = paste("Amino acid sequence:\n",
               substr(x$aaSeq, 1, 25), "...", substr(x$aaSeq, (nchar(x$aaSeq)-24), nchar(x$aaSeq)),
               "\n",
               sep = "")
    lines = c(lines, l5)
  }

  #addition if indel assessment
  if( "indel_likely" %in% names(x) ){
    if (x$indel_likely == TRUE){
      l6 = paste("The sequence likely contains an insertion or deletion.\nStop codon present: ",
                 x$stop_codons,
                 ", Amino acid PHMM score:",
                 x$aaScore,
                 sep = "")
    }else{
      l6 = paste("The sequence likely does not contain an insertion or deletion.\nStop codon present: ",
                 x$stop_codons,
                 ", Amino acid PHMM score:",
                 x$aaScore,
                 sep = "")
    }
    lines = c(lines, l6)
  }

  cat(lines,sep = "")
}
