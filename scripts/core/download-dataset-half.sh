#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-GT0-H.zip
wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-data-H.zip

unzip MiddEval3-GT0-H.zip -d ../../dataset
unzip MiddEval3-data-H.zip -d ../../dataset

rm MiddEval3-GT0-H.zip
rm MiddEval3-data-H.zip
