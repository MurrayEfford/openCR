###############################################################################
## package 'openCR'
## onAttach.R
## last changed 2018-03-11
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('openCR'), .openCRstuff$packageType)
    packageStartupMessage( "This is openCR ", version,
                           ". For overview type ?openCR and see openCR-vignette.pdf" )
}
