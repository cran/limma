read.imagene <- function(files,path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote="\"",...) {
#	Extracts an RG list from a series of Imagene analysis output files.
#	Imagene requires special treatment because red and green channel
#	intensities are in different files.
#	Gordon Smyth
#	14 Aug 2003.  Last modified 30 Dec 2006.

	if(is.null(dim(files))) {
		if(length(files)%%2==0)
			files <- matrix(files, ncol=2, byrow=TRUE)
		else
			stop("Odd number of files: should be two data files for each array")
	} else {
		files <- as.matrix(files)
		if(ncol(files) != 2) stop("Need a two column matrix of file names")
	}
	if(!is.null(ext)) files <- array(paste(files,ext,sep="."),dim(files))
	narrays <- nrow(files)
	if(is.null(names)) names <- paste(removeExt(files[,1]),removeExt(files[,2]),sep=".")

#	Read header information from first file to get nspots
	fullname <- files[1,1]
	if(!is.null(path)) fullname <- file.path(path,fullname)
	headers <- readImaGeneHeader(fullname)
	if(verbose) cat("Read header information\n")
	skip <- headers$NHeaderRecords
	FD <- headers$"Field Dimensions"
	if(is.null(FD)) stop("Can't find Field Dimensions in ImaGene header")
	nspots <- sum(apply(FD,1,prod))

#	Set signal and background estimation method
	if(is.null(columns)) {
		SM <- headers$"Measurement parameters"$"Segmentation Method"
		if(is.null(SM)) SM <- "none"
		if(SM=="auto")
			columns <- list(f="Signal Mean",b="Background Mean")
		else
			columns <- list(f="Signal Mean",b="Background Median")
	}
	if(is.null(columns$f) || is.null(columns$b)) stop("'columns' should have components 'f' and 'b'")

#	Initialize data object
	Y <- matrix(0,nspots,narrays)
	colnames(Y) <- names
	RG <- list(R=Y,G=Y,Rb=Y,Gb=Y,source="imagene",Field.Dimensions=FD)
	if(!is.null(wt.fun)) RG$weights <- Y

#	Other columns
	if(!is.null(other.columns)) {
		other.columns <- as.character(other.columns)
		if(length(other.columns)) {
			RG$other <- list()
			G.other.columns <- paste("G",other.columns)
			colnames(Y) <- removeExt(files[,1])
			for (a in G.other.columns) RG$other[[a]] <- Y 
			R.other.columns <- paste("R",other.columns)
			colnames(Y) <- removeExt(files[,2])
			for (a in R.other.columns) RG$other[[a]] <- Y 
		} else {
			other.columns <- NULL
		}
	}

	printer <- list(ngrid.r=FD[1,"Metarows"],ngrid.c=FD[1,"Metacols"],nspot.r=FD[1,"Rows"],nspot.c=FD[1,"Cols"])
	if(nrow(FD)==1) {
		RG$printer <- printer
	} else {
		printer$ngrid.r <- sum(FD[,"Metarows"])
		if(	all(printer$ngrid.c==FD[,"Metacols"]) &&
			all(printer$nspot.r==FD[,"Rows"]) &&
			all(printer$nspot.c==FD[,"Cols"]) ) RG$printer <- printer
	}

#	Now read data
	for (i in 1:narrays) {

#		Green channel
		fullname <- files[i,1]
		if(!is.null(path)) fullname <- file.path(path,fullname)
		if(i > 1) {
			headers <- readImaGeneHeader(fullname)
			if(any(headers[["Field Dimensions"]] != RG$Field.Dimensions))
				stop(paste("Field dimensions of array",i,"not same as those of first array"))
			skip <- headers$NHeaderRecords
		}
		obj<- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,check.names=FALSE,comment.char="",fill=TRUE,nrows=nspots,...)
		if(verbose) cat(paste("Read",fullname,"\n"))
		if(i==1) RG$genes <- obj[,c("Field","Meta Row","Meta Column","Row","Column","Gene ID")]
		RG$G[,i] <- obj[,columns$f]
		RG$Gb[,i] <- obj[,columns$b]
		if(!is.null(other.columns)) {
			for (j in 1:length(other.columns)) RG$other[[G.other.columns[j]]][,i] <- obj[,other.columns[j]]
		}

#		Red channel		
		fullname <- files[i,2]
		if(!is.null(path)) fullname <- file.path(path,fullname)
		obj<- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,check.names=FALSE,comment.char="",fill=TRUE,nrows=nspots,...)
		if(verbose) cat(paste("Read",fullname,"\n"))
		RG$R[,i] <- obj[,columns$f]
		RG$Rb[,i] <- obj[,columns$b]
		if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
		if(!is.null(other.columns)) {
			for (j in 1:length(other.columns)) RG$other[[R.other.columns[j]]][,i] <- obj[,other.columns[j]]
		}
	}
	new("RGList",RG)
}
