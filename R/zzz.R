#  ZZZ.R

#.First.lib <- function(libname, pkgname, where)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 23 Oct 2004.
#{
#	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) {
#		winMenuAddItem("Vignettes","limma","limmaUsersGuide()")
#	}
#}

.onAttach <- function(libname, pkgname)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	23 Jan 2004. Last revised 23 Oct 2004.
{
	if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) winMenuAddItem("Vignettes","limma","limmaUsersGuide()")
}
