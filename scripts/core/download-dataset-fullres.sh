#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-GT0-F.zip
wget http://vision.middlebury.edu/stereo/submit3/zip/MiddEval3-data-F.zip

unzip MiddEval3-GT0-F.zip -d ../../dataset
unzip MiddEval3-data-F.zip -d ../../dataset

rm MiddEval3-GT0-F.zip
rm MiddEval3-data-F.zip
