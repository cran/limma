#  ZZZ.R

limmaguideURL <- function() paste("file://",system.file("doc","usersguide.html",package="limma"),sep="")

.First.lib <- function(libname, pkgname, where)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 20 Sep 2004.
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) {
		winMenuAddItem("Vignettes","limma","browseURL(limmaguideURL())")
	}
}

