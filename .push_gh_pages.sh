#!/bin/bash

# Modified from:
# https://rmflight.github.io/posts/2014/11/travis_ci_gh_pages.html

# Completely overwrites gh-pages
# would be also possible to pull and modify instead
rm -rf out; || exit 0;
mkdir out;

GH_REPO="@github.com/antagomir/netresponse.git"

FULL_REPO="https://$GH_TOKEN$GH_REPO"

for files in '*.tar.gz'; do
        tar xfz $files
done

cd out
git init
git config user.name "netresponse-travis"
git config user.email "travis"
#cp ../microbiome/inst/doc/vignette.html index.html
#cp ../microbiome/vignettes/vignette.html index.html
touch index.html

git add .
git commit -m "Deployed to github pages"
git push --force --quiet $FULL_REPO master:gh-pages
