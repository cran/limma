#  ZZZ.R

limmaUsersGuide <- function(view=TRUE)
#	Find and optionally view limma User's Guide
#	Gordon Smyth
#	25 Oct 2004.
{
	f <- system.file("doc","usersguide.pdf",package="limma")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}

.First.lib <- function(libname, pkgname, where)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 23 Oct 2004.
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) {
		winMenuAddItem("Vignettes","limma","limmaUsersGuide()")
	}
}

