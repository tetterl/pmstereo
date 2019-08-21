#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

QUARTERIMGS_DIR=../../dataset/MiddEval3/trainingQ
SCALESERIES_DIR=../../dataset/TeddyScales
TEDDY_DIR=$QUARTERIMGS_DIR/Teddy
# make sure the quarter size images exist
if [ -d $QUARTERIMGS_DIR ]; then
  mkdir -p $SCALESERIES_DIR
  for SCALE in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100
  do
    echo "Downscaling to $SCALE %"
    convert -resize $SCALE%  $TEDDY_DIR/im0.png $SCALESERIES_DIR/im0-$SCALE.png
    convert -resize $SCALE%  $TEDDY_DIR/im1.png $SCALESERIES_DIR/im1-$SCALE.png
  done
else
  echo "Did not find quarter size images. Please run download-dataset.sh first!"
fi
