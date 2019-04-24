#library(devtools)
#library(roxygen2)
#library(Rcpp)
#document(".")
#install(".")

# I think roxygen2-based building generates Rcpp functions, therefore
Rcpp::compileAttributes()

# end then outside R CMD INSTALL --no-multiarch --with-keep.source .
