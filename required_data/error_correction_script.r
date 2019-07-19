# Below is the script I used in dapr testing for reading in data,
# translating it, passing it through the two PHMMs and setting reading frame.
# pull out the relevant code bits and use them here to interface with an R data
# strucutre named: COI5P

#need the structure to be something like an s3 class and to store the following :
# $rawSeq	$ntBin	$inFrame	$aaSeq	$aaBin	$aaScore	$indelLikely


#As you do this, move it over to the R/ folder and document.

#########################################################################

#benchmark start
start_time = Sys.time()
print(paste("Start time:", start_time))


print('loading packages & functions')
#required packages
#install.packages('ape')
library(ape)
#install.packages('aphid')
library(aphid)
#install.packages('seqinr')
library(seqinr)
#install.packages('tidyverse')
library(tidyverse)
# for mclapply
library(parallel)

source('../R/deploy_PHMMs.r')
source('../R/translation.r')

#When you get to a situation where there are multiple PHMMs to choose from,
#have it set up so that there is a lablelled list that determines which model set to load

#this is the main model block, would be run 
#file_list = c('Elasmobranchii_test.txt','all_bold_coi_v2.tsv')
file_list = c('all_bold_coi.tsv')


nt_phmm_file = 'COI5P_nt.PHMM'
aa_phmm_file = 'COI5P_aa.PHMM'

test_front = '../data/test/introduced_errors_test_'
output_front = '../data/test/model_output_test_'

#for multithreading of the processes
#numCores = detectCores()

#get the list of test files we are trying to correct with the models
#its good to keep them separate because it will reveal any bias towards
#different taxonomic groups
test_files = c() 
outfiles = c()

for(f in file_list){
	test_in = paste(test_front, f, sep = '')
	test_files = c(test_files, test_in) 
	
	outname = paste(output_front, f, sep = '')
	outfiles = c(outfiles, outname)

	}

names(outfiles) = test_files

correct_data = function(input){
	print("FOR DEVELOPMENT:")
	print("making a copy of the original nt_error column, as some will be altered and others will not")
	input$nt_error_unchanged = input$nt_error

#	print("working on column: nt_error. This will need to be interfaced with proper input column for real data")
#	#convert to a DNAbin
#	input$seq_DNAbins = mclapply(input$nt_error, function(x){
#		individual_DNAbin(x)
#		}, mc.cores = numCores)

	print("working on column: nt_error. This will need to be interfaced with proper input column for real data")
	#convert to a DNAbin
	input$seq_DNAbins = lapply(input$nt_error, function(x){
		individual_DNAbin(x)
		})

#	print("running Viterbi algorithm for all nucleotide sequences")
#	input$ntViterbi_out = mclapply(1:length(input$seq_DNAbins), function(i){
#		Viterbi(nt_PHMM, input$seq_DNAbins[[i]], odds = FALSE)
#		}, mc.cores = numCores)

	print("running slower Viterbi version with no mclapply")
	input$ntViterbi_out = lapply(1:length(input$seq_DNAbins), function(i){
		Viterbi(nt_PHMM, input$seq_DNAbins[[i]], odds = FALSE)
		})

	#run this with trycatch to see the problems:
#	print('checking for sequences with large numbers of leading non-profile base pairs')
#	input$leading_inserts = mclapply(1:length(input$ntViterbi_out), function(i){
#		leading_ins(input$ntViterbi_out[[i]][['path']])
#		}, mc.cores = numCores)

	#run this with trycatch to see the problems:
	print('checking for sequences with large numbers of leading non-profile base pairs')
	input$leading_inserts = unlist(lapply(1:length(input$ntViterbi_out), function(i){
		leading_ins(input$ntViterbi_out[[i]][['path']])
		}))

	#this chunck needs to be tested on small scale
	print('re-running ntPHMM for sequences with large numbers of inserts')
	print('This overwrites existing DNAbin and Viterbi output with the information from the shortened sequence')
	for(i in which(input$leading_inserts == TRUE)){
		print(paste("Long seqeuence at index:", i))
		tryCatch({
			#I'm using a different syntax here to draw attention to the fact I'm overwrtiting
			#previously built data for these special cases
			input[['nt_error']][[i]] = set_frame(input$nt_error[[i]], input$ntViterbi_out[[i]][['path']])
			input[['seq_DNAbins']][[i]] = individual_DNAbin(input[['nt_error']][[i]])
			input[['ntViterbi_out']][[i]] = Viterbi(nt_PHMM, input$seq_DNAbins[[i]], odds = FALSE)
		}, error = function(msg){
			message(paste("Error for row number:", i))
			return(NA)
		})
	}

	#commented out during test
		#this chunck needs to be tested on small scale
	#	print('re-running ntPHMM for sequences with large numbers of inserts')
	#	print('This overwrites existing DNAbin and Viterbi output with the information from the shortened sequence')
	#	for(i in which(input$leading_inserts == TRUE)){
	#		print(paste("Long seqeuence at index:", i))
	#		#I'm using a different syntax here to draw attention to the fact I'm overwrtiting
	#		#previously built data for these special cases
	#		input[['nt_error']][[i]] = set_frame(input$nt_error[[i]], input$ntViterbi_out[[i]][['path']])
	#		input[['seq_DNAbins']][[i]] = individual_DNAbin(input[['nt_error']][[i]])
	#		input[['ntViterbi_out']][[i]] = Viterbi(nt_PHMM, input$seq_DNAbins[[i]], odds = FALSE)
	#	}

	#previous version - tested but not robust to instances with large number of leading bp
	#print("setting reading frame for sequences")
	#print("working on column: nt_error. This will need to be interfaced with proper column for real data")
#	print("setting reading frame for sequences")
#	print("working on column: nt_error. This will need to be interfaced with proper column for real data")
#	input$nt_inFrame = unlist(mclapply(1:length(input$nt_error), function(i){
#		set_frame(input$nt_error[[i]], input$ntViterbi_out[[i]][['path']])
#		}, mc.cores = numCores))

	print("setting reading frame for sequences")
	print("working on column: nt_error. This will need to be interfaced with proper column for real data")
	input$nt_inFrame = unlist(lapply(1:length(input$nt_error), function(i){
		set_frame(input$nt_error[[i]], input$ntViterbi_out[[i]][['path']])
		}))

	#save the nt log odds score to a column in the dataframe
#	print("recording the Viterbi logl scores")
#	input$ntPHMM_score = unlist(mclapply(input$ntViterbi_out, function(x){
#		x[['score']]
#		}, mc.cores = numCores))	

	#save the nt log odds score to a column in the dataframe
	print("recording the Viterbi logl scores")
	input$ntPHMM_score = unlist(lapply(input$ntViterbi_out, function(x){
		x[['score']]
		}))	

#	print("censored translation of the sequences in reading frame")
#	input$translated_AA = unlist(mclapply(input$nt_inFrame, function(x){
#		censored_translation(x, reading_frame = 1)
#		}, mc.cores = numCores))

	print("censored translation of the sequences in reading frame")
	input$translated_AA = unlist(lapply(input$nt_inFrame, function(x){
		censored_translation(x, reading_frame = 1)
		}))

	#could have an if clause initiate the following if the user passes a known aa trans table
#	print(paste("translation of the sequences in reading frame with known amino acid table:", trans_table))
#	input$translated_AA = lapply(input$nt_inFrame, function(x){
#		trans_DNA_str(x, frame = 1,  numcode = 5, NAstring = '-')
#		})

#	print('converting amino acid sequences to binary format')
#	input$seq_AAbins = mclapply(input$translated_AA, function(x){
#		individual_AAbin(x)
#		}, mc.cores = numCores)

	print('converting amino acid sequences to binary format')
	input$seq_AAbins = lapply(input$translated_AA, function(x){
		individual_AAbin(x)
		})

#	print("running Viterbi algorithm for all amino acid sequences")
#	input$aaViterbi_out = mclapply(input$seq_AAbins, function(x){
#		Viterbi(aa_PHMM, x, odds = FALSE)
#		}, mc.cores = numCores)

	print("running Viterbi algorithm for all amino acid sequences")
	input$aaViterbi_out = lapply(input$seq_AAbins, function(x){
		Viterbi(aa_PHMM, x, odds = FALSE)
		})

	#save the nt log odds score to a column in the dataframe
#	print("recording the Viterbi logl scores")
#	input$aaPHMM_score = unlist(mclapply(input$aaViterbi_out, function(x){
#		x[['score']]
#		}, mc.cores = numCores))	

	#save the nt log odds score to a column in the dataframe
	print("recording the Viterbi logl scores")
	input$aaPHMM_score = unlist(lapply(input$aaViterbi_out, function(x){
		x[['score']]
		}))	

	print("-----------------------")
	print("Below only for development and optimization purposes, save the paths for the nt and amino")
	print("use these to assess performance by identifying the correction positions and comparing to the inputs")
	print("final version should instead output the logl, the seq in reading frame, and optionally the adjusted sequences.")
	print("-----------------------")

	#get the paths:
#	input$ntPath = unlist(mclapply(1:length(input$ntViterbi_out), function(i){
#		paste(input$ntViterbi_out[[i]][['path']], collapse = '')
#		}, mc.cores = numCores))

	#get the paths:
	input$ntPath = unlist(lapply(1:length(input$ntViterbi_out), function(i){
		paste(input$ntViterbi_out[[i]][['path']], collapse = '')
		}))

#	input$aaPath = unlist(mclapply(1:length(input$ntViterbi_out), function(i){
#		paste(input$aaViterbi_out[[i]][['path']], collapse='')
#		}, mc.cores = numCores))

	input$aaPath = unlist(lapply(1:length(input$ntViterbi_out), function(i){
		paste(input$aaViterbi_out[[i]][['path']], collapse='')
		}))

	#get the adjusted sequences
	input$nt_adjusted_sequences = unlist(lapply(1:length(input$nt_error), function(i){
		adj_seq(input$nt_error[[i]], input$ntViterbi_out[[i]][['path']])
		}))

	return(input)
}


#f = test_files[1]
for(f in test_files){

	print("loading the models")
	#load the PHMMs that are to be used:
	#nucleotide
	nt_PHMM = readPHMM(nt_phmm_file)
	#aa
	aa_PHMM = readPHMM(aa_phmm_file)

	print(paste("reading in file:", f))
	all_data = read_tsv(f)

#	print("first run on only 1% of test data")
#	print("remove this when it is working to scale up")
#	dim(all_data)
#	subsample = sample(1:length(all_data$nt_error), size = length(all_data$nt_error)/100)
#	all_data = all_data[subsample,]
#	dim(all_data)
#	all_data = head(all_data)

	subsets = seq(from = 100, to=length(all_data[[1]]), by=100)
	subsets = rep(seq_len(ceiling(length(all_data[[1]]) / 300)), each = 300, length.out = length(all_data[[1]]))

	print("spitting large dataframe into batches of 300 rows")
	batches = split(all_data, f = subsets)

	#free up memory
	all_data = NULL

	for(i in unique(subsets)){
		print(paste("on batch", i, "of", length(unique(subsets))))
		corrected = correct_data(batches[[i]])		

		dropcols = c("seq_DNAbins","ntViterbi_out","seq_AAbins","aaViterbi_out")

		corrected = corrected[ , !(names(corrected) %in% dropcols)]
		# if still failing use this to write at this point and append to the 
		# output file, header only for first iteration:
		if (i == 1){
			print("writing to output file")
			write.table(corrected, file = outfiles[[f]],
					quote = FALSE, row.names = FALSE,  sep = "\t", append=FALSE, col.names = TRUE) 	

		}else{
			print("appending to output file")
			write.table(corrected, file = outfiles[[f]], 
					quote = FALSE, row.names = FALSE, sep = "\t", append=TRUE, col.names = FALSE) 	
		}
		#
		#after writing then set the given df to an int to free up memory
		#cant set to null as this reindexes everything

		batches[[i]] = i
	}

}



end_time = Sys.time()
time_dif = difftime( end_time, start_time, units = "secs")
print(paste("Start time:", start_time))
print(paste("End time:", end_time))
print(paste("Correction of sequences took:", as.integer(time_dif), "seconds"))
