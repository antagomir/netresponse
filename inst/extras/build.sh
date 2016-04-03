#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check --as-cran netresponse_1.3.16.tar.gz
/usr/bin/R CMD INSTALL netresponse_1.3.16.tar.gz
#/usr/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
