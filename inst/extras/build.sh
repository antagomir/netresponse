#~/bin/R-patched/bin/R CMD BATCH document.R4.2.1
~/bin/R-4.2.1/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-4.2.1/bin/R CMD check netresponse_1.57.1.tar.gz #--no-build-vignettes --no-examples
#~/bin/R-4.2.1/bin/R CMD check --as-cran netresponse_1.57.1.tar.gz
~/bin/R-4.2.1/bin/R CMD BiocCheck netresponse_1.57.1.tar.gz
~/bin/R-4.2.1/bin/R CMD INSTALL netresponse_1.57.1.tar.gz 
#~/bin/R-4.2.1/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
