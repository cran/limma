limmaguideURL <- function() paste("file://",.find.package("limma"),"/doc/usersguide.html",sep="")

.First.lib <- function(libname, pkgname, where)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 2 July 2004.
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) {
#		if(!length(grep("^Biobase$", .packages()))) winMenuAdd("Vignettes")
		winMenuAddItem("Vignettes", "limma", "browseURL(limmaguideURL())")
	}
}
