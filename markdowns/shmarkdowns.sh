#!/bin/bash


#R -e "rmarkdown::render('0_check_100snps.Rmd',output_file='output.pdf')"

for filename in *.Rmd; do
  R -e "rmarkdown::render('$filename',output_file='$filename.pdf')" &
done
