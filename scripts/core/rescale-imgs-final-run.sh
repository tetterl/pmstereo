#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

FULLIMGS_DIR=../../dataset/MiddEval3/trainingF
SCALESERIES_DIR=../../dataset/size-dependency-run
TEDDY_DIR=$FULLIMGS_DIR/Teddy
# make sure the quarter size images exist
if [ -d $FULLIMGS_DIR ]; then
  mkdir -p $SCALESERIES_DIR
  for SCALE in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100
  do
    echo "Downscaling to $SCALE %"
    convert -resize $SCALE%  $TEDDY_DIR/im0.png $SCALESERIES_DIR/im0-$SCALE.png
    convert -resize $SCALE%  $TEDDY_DIR/im1.png $SCALESERIES_DIR/im1-$SCALE.png
  done
else
  echo "Did not find full size images. Please run download-dataset-fullres.sh first!"
fi
