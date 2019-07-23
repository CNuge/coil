# building the data and functions I've created for manipulating CPOI-5P sequences
# into a generic s3 function


########################
# Initialization of the class


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


}



#helper function - this is the user facing part of the class
coi5p = function(str = character()){
	
	#coerce the input into a lower case character string
	# vector can be another acceptable input

	#if vector, paste them together 
	validate_coi5p(new_coi5p(str))

}




###########################
# Generics and methods

#' this is where we would document this function in detail
#' this is a method dispatch, its job is to find the specific 
#' implementation for this class
frame = function(x){
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

	if(leading_ins(coi_obj$ntPHMM_out[['path']])){
		trim_temp  = set_frame(coi_obj$raw, ntPHMMout[['path']])
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
translate = function(x){
  UseMethod("translate")
}


translate.coi5p = function(coi_obj, frame = 0,  trans_table = 0){
	if(trans_table == 0){
		coi_obj$aaSeq = censored_translation(coi_obj$framed))
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
indel_check = function(x){
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