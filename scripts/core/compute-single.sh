#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.
echo $3
if [ "$#" -lt 2 ]
  then
    echo "Illegal number of arguments supplied. Usage:"
    echo "compute-single.sh IMAGE_NAME RUN_RESULT_DIR"
    echo "or"
    echo "compute-single.sh IMAGE_NAME RUN_RESULT_DIR ADDITIONAL_ARGS"
    echo "Note that IMAGE_NAME must contain the partial path prefix {trainingQ, trainingF}"
    exit 1
fi
# Path to the patch-match binary
BIN=../../build/pm
#Directory with the im0, im1 ground truth images
GT_DIR=../../dataset/MiddEval3/$1
# Path to the directory where the patch-match binary put it's output
RUN_RESULT=$2

LEFT_IMG=$GT_DIR/im0.png
RIGHT_IMG=$GT_DIR/im1.png

#Additional arguments (all after the first 2 mandatory ones)
ADDITIONAL_ARGS=${@:3}

# Creates output directory if it does not exist
mkdir -p $RUN_RESULT

#Running pm without additional args
if [ "$#" -eq 2 ]
then
  # Store current git branch and commit info
  touch $RUN_RESULT/run-params.txt
  echo "Git commit hash:" > $RUN_RESULT/run-params.txt
  echo "$(git rev-parse HEAD)" >> $RUN_RESULT/run-params.txt 
  echo "Git branch:" >> $RUN_RESULT/run-params.txt
  echo "$(git branch | grep \* | cut -d ' ' -f2)" >>  $RUN_RESULT/run-params.txt

  #Run Patchmatch
  $BIN --left $LEFT_IMG --right $RIGHT_IMG --output $RUN_RESULT
  
fi
#Running pm with additional arguments
if [ "$#" -gt 2 ]
then
  # Store current git branch and commit info
  # As well as any additional arguments we passed to patchmatch
  touch $RUN_RESULT/run-params.txt
  echo "Git commit hash:" > $RUN_RESULT/run-params.txt
  echo "$(git rev-parse HEAD)" >> $RUN_RESULT/run-params.txt 
  echo "Git branch:" > $RUN_RESULT/run-params.txt
  echo "$(git branch | grep \* | cut -d ' ' -f2)" >>  $RUN_RESULT/run-params.txt
  echo "Additional args:" >> $RUN_RESULT/run-params.txt 
  echo $ADDITIONAL_ARGS >> $RUN_RESULT/run-params.txt

  #Run Patchmatch
  $BIN --left $LEFT_IMG --right $RIGHT_IMG --output $RUN_RESULT $ADDITIONAL_ARGS
fi

