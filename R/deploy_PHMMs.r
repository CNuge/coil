# No imports - these should be in the workspace
#library(ape)
#library(aphid)
#library(seqinr)
#library(tidyverse)


#turn a column of DNA sequence strings into individual DNAbin objects
individual_DNAbin = function(dna_string){
	return(as.DNAbin(strsplit(gsub('-', 'n', as.character(tolower(dna_string))),"")))
}

#turn a column of AA sequence strings into individual AAbin objects
individual_AAbin = function(aa_string){
	return(as.AAbin(strsplit(as.character(aa_string),"")))
}

#TODO - think about combining leading_ins and ins_front_trim into a single function
#at present they're highly similar

#a check for a large number of leading inserted bases,
#if this is the case, TRUE is returned and the PHMM
#should be run a second time on the truncated data.
leading_ins = function(seq_path){
	#iterate along the sequence, if we hit 5 insertions
	#before 5 matches then there is a problem
	matches = 0
	ins = 0
	for(x in seq_path){
		if (x == 1){
			matches = matches + 1
		} else if (x == 2){
			ins = ins + 1
		}

		if (matches == 5){
			return(FALSE)
		}
		if (ins == 5){
			return(TRUE)
		}
	}
}

# check sequence for an early large string of deletions, if it exists then
# return the starting index by which to slice the path and the string
ins_front_trim = function(path_out, search_scope = 15){

	if(sum(path_out[1:search_scope] == 2)>2){
		run_pos = min(which(path_out == 2))
		for(i in run_pos:length(path_out)){
			if(path_out[i] != 2){
				return(i)
			}
		}
	}

	return(FALSE)
}


#take an input sequence and get it into the reading frame for proper translation based on the
#nucleotide PHMM output. This is more conservative, just what is required to get the aaPHMM
#running, as opposed to making corrections to the nucleotides first
set_frame = function(org_seq , path_out){
	org_seq_vec = strsplit(tolower(org_seq), split='')[[1]]

	adj_for_dels = ins_front_trim(path_out)

	#If there are a large number (>2) of deletes in the first 15bp
	#then adjust the sequence so the end of the delete run is the starting point
	if(adj_for_dels != FALSE){
		org_seq_vec = org_seq_vec[adj_for_dels:length(org_seq_vec)]
		path_out = path_out[adj_for_dels:length(path_out)]
	}

	front = c()
	org_seq_start = 1
	org_seq_end = length(org_seq_vec)


	for( i in 1:length(path_out) ){
		#0 = D
		if( path_out[i] == 0 ){
			#there was a bp missing in the original seq
			#add a placeholder to the new seq, don't advance on original pointer
			front = c(front, '-')
		#2 = I
		}else if( path_out[i] == 2 ){
			#there was an extra bp at the front of the sequence
			#not represented in the PHMMs, skip this bp in the original seq
			org_seq_start = org_seq_start + 1
		#for the first match seen, check to make sure it isn't a single match dangling
		#in a sea of inserts or deletes
		}else if( path_out[i] == 1 ){
			if (path_out[i+1] == 2 & path_out[i+2] == 2){
				org_seq_start = org_seq_start + 1
			}else if (path_out[i+1] == 0 & path_out[i+2] == 0){
				front = c(front, '-')
			}else{
				break
			}
		}
	}

	for( i in length(path_out):1 ){
		#0 = D
		if( path_out[i] == 0 ){
			#there was a bp missing at the back of the original seq
			#just continue to the next position
			next
		#2 = I
		}else if( path_out[i] == 2 ){
			#there is an extra base pair at the end of the original sequence
			#not represented in the PHMMs, remove this trailing end of the original sequence
			org_seq_end = org_seq_end - 1
		}else if ( path_out[i] == 1){
			#in either of these instances we want to trim the last bp, as its
			#dangling in a sea of inserts or deletes and likely a random profile match
			if (path_out[i-1] == 2 & path_out[i-2] == 2){
				org_seq_end = org_seq_end - 1
			}else if (path_out[i-1] == 0 & path_out[i-2] == 0){
				org_seq_end = org_seq_end - 1
			}else{
				break
			}
		}
	}

	return(paste(c(front,org_seq_vec[org_seq_start:org_seq_end]),collapse= ""))
}




