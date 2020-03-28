#~/bin/ CMD BATCH document.R
~/bin/R-3.6.2/bin/R CMD build ../../
~/bin/R-3.6.2/bin/R CMD check --as-cran netresponse_1.47.2.tar.gz
~/bin/R-3.6.2/bin/R CMD INSTALL netresponse_1.47.2.tar.gz
#~/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
