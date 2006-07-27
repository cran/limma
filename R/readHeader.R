#  readHeader.R
#  Read header information from top of file
#  for various microarray image analysis formats

readGenericHeader <- function(file,columns,sep="\t")
#	Find line on which numeric data starts in microarray data file
#	using only names of data columns
#	Gordon Smyth
#	10 March 2006. Last modified 5 April 2006.
{
	if(missing(columns) || !length(columns)) stop("must specify column headings to find")
	columns <- protectMetachar(as.character(columns))
	if(!length(columns)) stop("column headings must be specified")
	con <- file(file, "r")	
 	on.exit(close(con))
 	out <- list()
 	Found <- FALSE
 	i <- 0
 	repeat {
 		i <- i+1
		txt <- readLines(con,n=1)
		if(!length(txt)) stop("Specified column headings not found in file")
		Found <- TRUE
		for(a in columns) Found <- Found && length(grep(a,txt))
		if(Found) break
	}
	out$NHeaderRecords <- i-1
	out$ColumnNames <- strsplit(txt,split=sep)[[1]]
	out
}

readGPRHeader <- function(file) {
#	Extracts header information from a GenePix Results (GPR) file
#	Gordon Smyth
#	4 October 2003.  Last modified 25 June 2004.

	con <- file(file, "r")	
 	on.exit(close(con))
 	if(substring(readLines(con,n=1),1,3) != "ATF") stop("File is not in Axon Text File (ATF) format")
	nfields <- as.numeric(strsplit(readLines(con,n=1),split="\t")[[1]])
	out <- list()
	for (i in 1:nfields[1]) {
		txt <- strsplit(sub("\"$","",sub("^\"","",readLines(con,n=1))),split="=")[[1]]
		out[[i]] <- txt[2]
		names(out)[i] <- txt[1]
	}
	out$NHeaderRecords <- nfields[1]+2
	out$NDataFields <- nfields[2]
	out
}

readImaGeneHeader <- function(file)
#	Extracts header information from an Imagene analysis output file
#	Gordon Smyth
#	14 Aug 2003.  Last modified 27 July 2006.
{
	con <- file(file, "r")	
 	on.exit(close(con))
 	out <- list(Header=list())
 	iline <- 0
 	CompName <- character(0)
 	repeat {
 		txt <- trimWhiteSpace(readLines(con,n=1))
 		if(!length(txt)) break
		iline <- iline+1
		if (length(grep("^Begin",txt))) {
			CompName <- c(CompName,sub("^Begin ","",txt))
			if(CompName[1] != "Header") stop("Begin Header missing or mis-matched Begin/Ends")
			out[[CompName]] <- list()
		} else {
			if (length(grep("^End",txt))) {
				CompName <- CompName[-length(CompName)]
	 			if(txt=="End Header") {
	 				out <- out$Header
	 				out$NHeaderRecords <- iline+1
	 				if(!is.null(out[["Field Dimensions"]])) {
	 					cn <- out[[c("Field Dimensions","Field")]]
	 					rn <- names(out[["Field Dimensions"]])[-1]
	 					x <- matrix(as.integer(0),length(rn),length(cn),dimnames=list("Field"=rn,cn))
	 					for(a in rn) x[a,] <- as.integer(out[[c("Field Dimensions",a)]])
	 					out[["Field Dimensions"]] <- x
	 				}
	 				if(!is.null(out[["Alerts"]])) {
	 					cn <- out[[c("Alerts","Control Type")]]
	 					rn <- names(out[["Alerts"]])[-1]
	 					x <- matrix("",length(rn),length(cn),dimnames=list("Control Type"=rn,cn))
	 					for(a in rn) x[a,] <- out[[c("Alerts",a)]]
	 					out[["Alerts"]] <- x
	 				}
					return(out)
	 			}
			} else {
				split <- "\t"
				if(!length(grep(split,txt))) split <- ": "
 				txtsplit <- strsplit(txt,split=split)[[1]]
 				n <- length(txtsplit)
 				if(n>0) {
 					FieldName <- txtsplit[1]
 					if(n>1) FieldValue <- txtsplit[-1] else FieldValue <- ""
					a <- c(CompName,FieldName)					
					out[[a]] <- FieldValue
				}
			}
		}
 	}
 	warning("End of file encountered before End Header")
 	out
}

# readQuantArrayHeader <- function(file)
# #	Extracts header information from a QuantArray output file
# #	Gordon Smyth
# #	20 May 2006.
# {
# 	con <- file(file, "r")	
#  	on.exit(close(con))
#  	out <- list(Header=list())
#  	iline <- 0
#  	CompName <- character(0)
#  	repeat {
#  		txt <- readLines(con,n=1)
#  		if(!length(txt)) break
# 		iline <- iline+1
# 		if (length(grep("^Begin",txt))) {
# 			CompName <- c(CompName,sub("^Begin ","",txt))
# 			out[[CompName]] <- list()
# 		} else {
# 			if (length(grep("^End",txt))) {
# 				CompName <- CompName[-length(CompName)]
# 	 			if(txt=="End Header") {
# 	 				out <- out$Header
# 	 				out$NHeaderRecords <- iline+1
# 					return(out)
# 	 			}
# 			} else {
# 				split <- "\t"
# 				if(!length(grep(split,txt))) split <- ": "
#  				txtsplit <- strsplit(txt,split=split)[[1]]
#  				n <- length(txtsplit)
#  				if(n>0) {
#  					FieldName <- txtsplit[1]
#  					if(n>1) FieldValue <- txtsplit[-1] else FieldValue <- ""
# 					a <- c(CompName,FieldName)					
# 					out[[a]] <- FieldValue
# 				}
# 			}
# 		}
#  	}
#  	warning("End of file encountered before End Header")
#  	out
# }

readSMDHeader <- function(file)
#	Read header information from a Stanford Microarray Database (SMD) raw data file
#	Gordon Smyth
#	3 June 2004
{
	con <- file(file, "r")	
 	on.exit(close(con))
 	out <- list()
 	i <- 0
 	repeat {
 		txt <- readLines(con,n=1)
 		if(txt=="") stop("Not SMD data file: input stopped at blank line")
 		if(substr(txt,1,1)=="!") {
 	 		i <- i+1
 	 		txtsplit <- strsplit(substr(txt,2,500),split="=")[[1]]
			out[[i]] <- txtsplit[2]
			names(out)[i] <- txtsplit[1]
 		} else {
 			break
 		}
 	}
 	out$NHeaderRecords <- i
 	out
}

