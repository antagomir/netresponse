#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check netresponse_1.21.14.tar.gz
/usr/bin/R CMD INSTALL netresponse_1.21.14.tar.gz
#/usr/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
