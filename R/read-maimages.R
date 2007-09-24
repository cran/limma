#	READ IMAGE ANALYSIS FILES INTO RGList

read.maimages <- function(files=NULL,source="generic",path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,annotation=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
#	Extracts an RG list from a series of image analysis output files
#	Gordon Smyth. 
#	1 Nov 2002.  Last revised 3 Feb 2007.
{
#	Begin checking input arguments

	if(is.null(files)) {
		if(is.null(ext))
			stop("Must specify input files")
		else {
			extregex <- paste("\\.",ext,"$",sep="")
			files <- dir(path=ifelse(is.null(path),".",path),pattern=extregex)
			files <- sub(extregex,"",files)
		}
	}

	source <- match.arg(source,c("generic","agilent","arrayvision","bluefuse","genepix","genepix.median","genepix.custom","imagene","quantarray","scanarrayexpress","smd.old","smd","spot","spot.close.open"))
#	source2 is the source type with qualifications removed
	source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
	if(is.null(quote)) if(source=="agilent") quote <- "" else quote <- "\""
	if(source2=="imagene") return(read.imagene(files=files,path=path,ext=ext,names=names,columns=columns,other.columns=other.columns,wt.fun=wt.fun,verbose=verbose,sep=sep,quote=quote,...))

	if(is.data.frame(files)) {
		targets <- files
		files <- files$FileName
		if(is.null(files)) stop("targets frame doesn't contain FileName column")
		if(is.null(names)) names <- targets$Label
	} else {
		targets <- NULL
	}

	slides <- as.vector(as.character(files))
	if(!is.null(ext)) slides <- paste(slides,ext,sep=".")
	nslides <- length(slides)

	if(is.null(names)) names <- removeExt(files)

	if(is.null(columns)) {
		if(source2=="generic") stop("must specify columns for generic input")
		columns <- switch(source,
			agilent = list(G="gMeanSignal",Gb="gBGMedianSignal",R="rMeanSignal",Rb="rBGMedianSignal"),
			arrayvision = list(G="MTM Dens - Levels",Gb="Bkgd",R="MTM Dens - Levels",Rb="Bkgd"),
			bluefuse = list(G="AMPCH1",R="AMPCH2"),
			genepix = list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"),
			genepix.median = list(R="F635 Median",G="F532 Median",Rb="B635 Median",Gb="B532 Median"),
			genepix.custom = list(R="F635 Mean",G="F532 Mean",Rb="B635",Gb="B532"),
			quantarray = list(R="ch2 Intensity",G="ch1 Intensity",Rb="ch2 Background",Gb="ch1 Background"),
			scanarrayexpress = list(G="Ch1 Mean",Gb="Ch1 B Median",R="Ch2 Mean",Rb="Ch2 B Median"),
			smd.old = list(G="CH1I_MEAN",Gb="CH1B_MEDIAN",R="CH2I_MEAN",Rb="CH2B_MEDIAN"),
			smd = list(G="Ch1 Intensity (Mean)",Gb="Ch1 Background (Median)",R="Ch2 Intensity (Mean)",Rb="Ch2 Background (Median)"),
			spot = list(R="Rmean",G="Gmean",Rb="morphR",Gb="morphG"),
			spot.close.open = list(R="Rmean",G="Gmean",Rb="morphR.close.open",Gb="morphG.close.open"),
			NULL
		)
	} else {
		if(!is.list(columns)) stop("columns must be a list")
		names(columns)[names(columns)=="Gf"] <- "G"
		names(columns)[names(columns)=="Rf"] <- "R"
		if(is.null(columns$G) || is.null(columns$R)) stop("columns must specify foreground G and R")
		if(!all(names(columns) %in% c("G","R","Gb","Rb"))) warning("non-standard columns specified")
	}
	cnames <- names(columns)
	
	if(is.null(annotation)) annotation <- switch(source,
		agilent = c("Row","Col","Start","Sequence","SwissProt","GenBank","Primate","GenPept","ProbeUID","ControlType","ProbeName","GeneName","SystematicName","Description"),
		arrayvision = "Spot labels",
		bluefuse = c("ROW","COL","SUBGRIDROW","SUBGRIDCOL","BLOCK","NAME","ID"),   
		genepix=,genepix.median=,genepix.custom = c("Block","Row","Column","ID","Name"),
		quantarray= c("Array Row","Array Column","Row","Column","Name"),
		scanarrayexpress = c("Array Row","Array Column","Spot Row","Spot Column"), 	
		smd = c("Spot","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","Locuslink ID","Name","Sequence Type","X Grid Coordinate (within sector)","Y Grid Coordinate (within sector)","Sector","Failed","Plate Number","Plate Row","Plate Column","Clone Source","Is Verified","Is Contaminated","Luid"),
		smd.old = c("SPOT","NAME","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","SUID"),
		NULL
	)

#	End checking input arguments

#	Read first file to get nspots
	fullname <- slides[1]
	if(!is.null(path)) fullname <- file.path(path,fullname)
	required.col <- unique(c(annotation,unlist(columns),other.columns))
	text.to.search <- if(is.null(wt.fun)) "" else deparse(wt.fun)
	switch(source2,
		quantarray = {
		firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
		skip <- grep("Begin Data",firstfield)
		if(length(skip)==0) stop("Cannot find \"Begin Data\" in image output file")
		nspots <- grep("End Data",firstfield) - skip -2
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,nrows=nspots,flush=TRUE,...)
	}, arrayvision = {
		skip <- 1
		cn <- scan(fullname,what="",sep=sep,quote=quote,skip=1,nlines=1,quiet=TRUE,allowEscape=FALSE)
		fg <- grep(" Dens - ",cn)
		if(length(fg) != 2) stop(paste("Cannot find foreground columns in",fullname))
		bg <- grep("^Bkgd$",cn)
		if(length(bg) != 2) stop(paste("Cannot find background columns in",fullname))
#		Note that entries for columns for ArrayVision are now numeric
		columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,flush=TRUE,...)
#		obj <- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,as.is=TRUE,check.names=FALSE,fill=TRUE,comment.char="",flush=TRUE,...)
		fg <- grep(" Dens - ",names(obj))
		bg <- grep("^Bkgd$",names(obj))
		columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
		nspots <- nrow(obj)
	}, bluefuse = {
		skip <- readGenericHeader(fullname,columns=c(columns$G,columns$R))$NHeaderRecords
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,flush=TRUE,...)
		nspots <- nrow(obj)
	}, genepix = {
		h <- readGPRHeader(fullname)
		if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
		skip <- h$NHeaderRecords
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,flush=TRUE,...)
		nspots <- nrow(obj)
	}, smd = {
		skip <- readSMDHeader(fullname)$NHeaderRecords
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,flush=TRUE,...)
		nspots <- nrow(obj)
	}, {
		skip <- readGenericHeader(fullname,columns=columns,sep=sep)$NHeaderRecords
		obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,as.is=TRUE,fill=TRUE,flush=TRUE,...)
		nspots <- nrow(obj)
	})

#	Initialize RG list object (object.size for matrix of NAs is smaller)
	Y <- matrix(NA,nspots,nslides)
	colnames(Y) <- names
	RG <- columns
	for (a in cnames) RG[[a]] <- Y
	if(!is.null(wt.fun)) RG$weights <- Y
	if(is.data.frame(targets)) {
		RG$targets <- targets
	} else {
		RG$targets <- data.frame(FileName=I(files),row.names=names)
	}

#	Set annotation columns
	if(!is.null(annotation)) {
		j <- match(annotation,colnames(obj),0)
		if(any(j>0)) RG$genes <- data.frame(obj[,j,drop=FALSE],check.names=FALSE)
	}

	RG$source <- source

#	Set printer layout, if possible
	if(source2=="agilent") {
		if(!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
			nr <- length(unique(RG$genes$Row))
			nc <- length(unique(RG$genes$Col))
			if(nspots==nr*nc) RG$printer <- list(ngrid.r=1,ngrid.c=1,nspot.r=nr,nspot.c=nc)
		}
	}
	if(source2=="genepix") {
		if(!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && !is.null(RG$genes$Column)) {
			RG$printer <- getLayout(RG$genes,guessdups=FALSE)
			nblocks <- RG$printer$ngrid.r*RG$printer$ngrid.c
			if(!is.na(nblocks) && (nblocks>1) && !is.null(obj$X)) {
				blocksize <- RG$printer$nspot.r*RG$printer$nspot.c
				i <- (1:(nblocks-1))*blocksize
				ngrid.r <- sum(obj$X[i] > obj$X[i+1]) + 1
				if(!is.na(ngrid.r) && nblocks%%ngrid.r==0) {
					RG$printer$ngrid.r <- ngrid.r
					RG$printer$ngrid.c <- nblocks/ngrid.r
				} else {
					warning("Can't determine number of grid rows")
					RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
				}
			}
		}
	}

#	Other columns
	if(!is.null(other.columns)) {
		other.columns <- as.character(other.columns)
		j <- match(other.columns,colnames(obj),0)
		if(any(j>0)) {
			other.columns <- colnames(obj)[j]
			RG$other <- list()
			for (j in other.columns) RG$other[[j]] <- Y 
		} else {
			other.columns <- NULL
		}
	}

#	Read remainder of files
	for (i in 1:nslides) {
		if(i > 1) {
			fullname <- slides[i]
			if(!is.null(path)) fullname <- file.path(path,fullname)
			switch(source2,
			   quantarray = {
					firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
					skip <- grep("Begin Data", firstfield)
				},  arrayvision = {
					skip <- 1
				}, genepix = {
					skip <- readGPRHeader(fullname)$NHeaderRecords
				}, smd = {
					skip <- readSMDHeader(fullname)$NHeaderRecords
				}, {
					skip <- readGenericHeader(fullname,columns=columns)$NHeaderRecords
				})
			if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
			obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,as.is=TRUE,quote=quote,fill=TRUE,nrows=nspots,flush=TRUE,...)
		}
		for (a in cnames) RG[[a]][,i] <- obj[,columns[[a]]]
		if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
		if(!is.null(other.columns)) for (j in other.columns) {
			RG$other[[j]][,i] <- obj[,j] 
		}
		if(verbose) cat(paste("Read",fullname,"\n"))
	}
	new("RGList",RG)
}

read.columns <- function(file,required.col=NULL,text.to.search="",sep="\t",quote="\"",skip=0,fill=TRUE,blank.lines.skip=TRUE,comment.char="",allowEscapes=FALSE,...)
#	Read specified columns from a delimited text file with header line
#	Gordon Smyth
#	3 Feb 2007. Last modified 1 May 2007.
{
#	Default is to read all columns
	if(is.null(required.col)) return(read.table(file=file,header=TRUE,check.names=FALSE,sep=sep,quote=quote,skip=skip,fill=fill,blank.lines.skip=blank.lines.skip,comment.char=comment.char,allowEscapes=allowEscapes,...))

	text.to.search <- as.character(text.to.search)

#	Read header to get column names
    allcnames <- scan(file,what="",sep=sep,quote=quote,nlines=1,quiet=TRUE,skip=skip,strip.white=TRUE,blank.lines.skip=blank.lines.skip,comment.char=comment.char,allowEscapes=allowEscapes)
	ncn <- length(allcnames)
	colClasses <- rep("NULL",ncn)

#	Are required columns in the header?
	if(is.numeric(required.col)) {
		colClasses[required.col] <- NA
	} else {
		colClasses[allcnames %in% as.character(required.col)] <- NA
	}

#	Search for column names in text
	if(any(text.to.search != "")) for (i in 1:ncn) {
		if(length(grep(protectMetachar(allcnames[i]),text.to.search))) colClasses[i] <- NA
	}

#	Is there a leading column of row.names without a header?
    secondline <- scan(file,what="",sep=sep,quote=quote,nlines=1,quiet=TRUE,skip=skip+1,strip.white=TRUE,blank.lines.skip=blank.lines.skip,comment.char=comment.char,allowEscapes=allowEscapes)
	if(length(secondline) > ncn) colClasses <- c(NA,colClasses)

#	Read specified columns
	read.table(file=file,header=TRUE,col.names=allcnames,check.names=FALSE,colClasses=colClasses,sep=sep,quote=quote,skip=skip,fill=fill,blank.lines.skip=blank.lines.skip,comment.char=comment.char,allowEscapes=allowEscapes,...)
}

