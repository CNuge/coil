# building the data and functions I've created for manipulating CPOI-5P sequences
# into a generic s3 function



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