# no imports - this should be in the workspace
#library(seqinr)


#' Censored Translation
#' Translate a DNA sequence using the censored translation table,
#' this translates codons for which the amino acids is unambigious across the
#' animal kingdom, and does not translate those for which the amino acid varies
#' but rather outputs a ? in the string
#` Censored translation table:
#`            FFLLSSSSYY?*CCWWLLLLPPPPHHQQRRRRII?MTTTTNN?KSS??VVVVAAAADDEEGGGG
#`   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#`   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#`   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
translate_codon = function(codon){

	trans_table = list(
		'TTT' = 'F', 'TTC' = 'F', 'TTA' = 'L', 'TTG' = 'L',
		'TCT' = 'S', 'TCC' = 'S', 'TCA' = 'S', 'TCG' = 'S',
		'TAT' = 'Y', 'TAC' = 'Y', 'TAA' = '?', 'TAG' = '*',
		'TGT' = 'C', 'TGC' = 'C', 'TGA' = 'W', 'TGG' = 'W',
		'CTT' = 'L', 'CTC' = 'L', 'CTA' = 'L', 'CTG' = 'L',
		'CCT' = 'P', 'CCC' = 'P', 'CCA' = 'P', 'CCG' = 'P',
		'CAT' = 'H', 'CAC' = 'H', 'CAA' = 'Q', 'CAG' = 'Q',
		'CGT' = 'R', 'CGC' = 'R', 'CGA' = 'R', 'CGG' = 'R',
		'ATT' = 'I', 'ATC' = 'I', 'ATA' = '?', 'ATG' = 'M',
		'ACT' = 'T', 'ACC' = 'T', 'ACA' = 'T', 'ACG' = 'T',
		'AAT' = 'N', 'AAC' = 'N', 'AAA' = '?', 'AAG' = 'K',
		'AGT' = 'S', 'AGC' = 'S', 'AGA' = '?', 'AGG' = '?',
		'GTT' = 'V', 'GTC' = 'V', 'GTA' = 'V', 'GTG' = 'V',
		'GCT' = 'A', 'GCC' = 'A', 'GCA' = 'A', 'GCG' = 'A',
		'GAT' = 'D', 'GAC' = 'D', 'GAA' = 'E', 'GAG' = 'E',
		'GGT' = 'G', 'GGC' = 'G', 'GGA' = 'G', 'GGG' = 'G'
		)

	if(nchar(codon) != 3){
		return('-')
	}

	if(grepl('-', codon,  fixed = TRUE)){
		return('-')
	}

	if(grepl('N', toupper(codon),  fixed = TRUE)){
		return('-')
	}

	return(trans_table[[toupper(codon)]])
}


censored_translation = function(dna_str, reading_frame = 1){
	#reading frame = 1 means the first bp in the string is the start of the
	#first codon, can pass 1, 2 or 3. For 2 and 3 the first 1 and 2 bp will be
	#dropped from translation respectively
	num_bp = nchar(dna_str)

	codons = seq(reading_frame, num_bp, by=3)

	codon_vec = sapply(codons, function(x) {
						substr(dna_str, x, x+2)
						})

	aa_str = paste(lapply(codon_vec, translate_codon), collapse= "")

	return(aa_str)
}


#`
#`
#`
#`
#`
trans_dna = function(dna_str, frame = 0,  trans_table = 0){
	if(trans_table == 0){
		if(frame != 0){
		  return(censored_translation(substring(dna_str, frame+1)))
		}
	  return(censored_translation(dna_str))
	}else{
		#split the DNA string into a vector, all characters to lower case
		dna_list = strsplit(gsub('-', 'n', as.character(tolower(dna_str))),"")
		dna_vec = dna_list[[1]]
		#translate using the designated numcode, returns a vector of AAs
		aa_vec = seqinr::translate(dna_vec, frame = frame, numcode=trans_table, ambiguous= TRUE, NAstring = '-')

		aa_str = paste(aa_vec, collapse= "")
		return(aa_str)
	}
}


#path would be relative to the library design
#translation_table_data = read.table('../required_data/family_tanslation_table.tsv' ,
#								header = TRUE, sep = '\t', stringsAsFactors = FALSE)


#` determine the translation table to use for a given phylogenetic group
#` data stored down to family level.
#` relies on the above having been run so that the df is in the workspace and accessable
#'
#'@param x a taxonomic designation (family order  class  phylum) first letter capitilizaed
#'
#'@return an integer indicating the correct translation table (give bcbi link here)
#'@examples
#'
#'@export
which_trans_table = function(x) {
  trans_df$trans_table[trans_df$taxon == x]

}

#
