#!/bin/sh
Rscript -e 'bookdown::render_book("_output.yml", output_format = "all")'
rm -rf book-output
git clone -b gh-pages \
  https://github.com/LeonieHagitte/Thesis.git \
  book-output
cd book-output
git rm -rf *
cp -r ../_book/* ./
git add --all *
git commit -m "Update the thesis"
git push -q origin gh-pages
