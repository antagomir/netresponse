/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check netresponse_1.17.11.tar.gz
/usr/bin/R CMD INSTALL netresponse_1.17.11.tar.gz
#/usr/bin/R CMD BiocCheck netresponse_1.17.1.tar.gz 


