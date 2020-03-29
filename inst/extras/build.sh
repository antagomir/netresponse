#~/bin/R-patched/bin/R CMD BATCH document.R
~/bin/R-patched/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-patched/bin/R CMD check netresponse_1.47.3.tar.gz #--no-build-vignettes --no-examples
#~/bin/R-patched/bin/R CMD check --as-cran netresponse_1.47.3.tar.gz
~/bin/R-patched/bin/R CMD BiocCheck netresponse_1.47.3.tar.gz
~/bin/R-patched/bin/R CMD INSTALL netresponse_1.47.3.tar.gz 
#~/bin/R-patched/bin/R CMD BiocCheck netresponse_1.17.13.tar.gz 
