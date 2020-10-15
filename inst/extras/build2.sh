#~/bin/ CMD BATCH document.R
~/bin/R-4.0.0/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
# ~/bin/R-4.0.0/bin/R CMD check --as-cran netresponse_1.49.1.tar.gz
~/bin/R-4.0.0/bin/R CMD BiocCheck netresponse_1.49.1.tar.gz
~/bin/R-4.0.0/bin/R CMD INSTALL netresponse_1.49.1.tar.gz
#~/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
