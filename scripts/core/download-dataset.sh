#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-GT0-Q.zip
wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-data-Q.zip

unzip MiddEval3-GT0-Q.zip -d ../../dataset
unzip MiddEval3-data-Q.zip -d ../../dataset

rm MiddEval3-GT0-Q.zip
rm MiddEval3-data-Q.zip
