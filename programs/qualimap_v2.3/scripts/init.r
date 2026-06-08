
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

required.packages <- c('optparse', 'NOISeq', 'XML', 'Rsamtools', 'Repitools', 'rtracklayer' )

for (package in required.packages) {
  if (!is.installed(package)) {
    print (paste("ERROR! Package :", package, ": is missing!"))
  }
}

