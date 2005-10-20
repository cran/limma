read.imagene <- function(files,path=NULL,ext=NULL,names=NULL,columns=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote="\"",...) {
#	Extracts an RG list from a series of Imagene analysis output files.
#	Imagene requires special treatment because red and green channel
#	intensities are in different files.
#	Gordon Smyth
#	14 Aug 2003.  Last modified 18 September 2005.

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
	if(is.null(columns)) columns <- list(f="Signal Mean",b="Background Median")
	if(is.null(columns$f) || is.null(columns$b)) stop("'columns' should have components 'f' and 'b'")
	if(is.null(names)) names <- removeExt(files[,1])

#	Read header information from first file to get nspots
	fullname <- files[1,1]
	if(!is.null(path)) fullname <- file.path(path,fullname)
	headers <- readImaGeneHeader(fullname)
	if(verbose) cat("Read header information\n")
	skip <- headers$NHeaderRecords
	FD <- headers[["Field Dimensions"]]
	nspots <- prod(FD)

#	Now read data
	Y <- matrix(0,nspots,narrays)
	colnames(Y) <- names
	RG <- list(R=Y,G=Y,Rb=Y,Gb=Y,Image.Program="ImaGene",Field.Dimensions=FD)
	if(!is.null(wt.fun)) RG$weights <- Y
	if(nrow(FD)==1) {
		RG$printer <- list(ngrid.r=FD[1,"Metarows"],ngrid.c=FD[1,"Metacols"],nspot.r=FD[1,"Rows"],nspot.c=FD[1,"Cols"])
	}
	for (i in 1:narrays) {
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
		fullname <- files[i,2]
		if(!is.null(path)) fullname <- file.path(path,fullname)
		obj<- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,check.names=FALSE,comment.char="",fill=TRUE,nrows=nspots,...)
		if(verbose) cat(paste("Read",fullname,"\n"))
		RG$R[,i] <- obj[,columns$f]
		RG$Rb[,i] <- obj[,columns$b]
		if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
	}
	new("RGList",RG)
}
