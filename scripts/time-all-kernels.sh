#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.
#We require 3 parameters (if running all kernels)
#else more
if [ "$#" -lt 4 ]
  then
    echo "Illegal number of arguments supplied. Usage:"
    echo "./time-all-kernels.sh {trainingQ, trainingF} {image name} EVAL RUN_SUBSET"
    echo "EVAL = 1 -> run evaluation"
    echo "EVAL = 0 -> just measure timings"
    echo "RUN_SUBSET = 1 -> run subset of kernels specified"
    echo "RUN_SUBSET = 0 -> run all kernels"
    echo "example:"
    echo "./time-all-kernels.sh trainingQ Teddy 0 0 "
    exit 1
fi

IMAGE_DIR=$1 #choose trainingF for fullrs
IMAGE_NAME=$2 
RUN_RESULT=../runs/$IMAGE_DIR/$IMAGE_NAME

# Create output directory if it does not exist
mkdir -p $RUN_RESULT

#Eval paths and variables
EVAL_BIN=../build/eval

#Directory with the im0, im1 ground truth images
GT_DIR=../dataset/MiddEval3/$IMAGE_DIR/$IMAGE_NAME
LEFT_IMG=$GT_DIR/im0.png
RIGHT_IMG=$GT_DIR/im1.png
PFM_FILE=$RUN_RESULT/disp1.pfm #output pfm
GT_FILE=$GT_DIR/disp0GT.pfm #gt input pfm
MASK_FILE=$GT_DIR/mask0nocc.png

#If this is set to 0, ALL kernels are run!
RUN_SUBSET_OF_KERNELS=$4 #{True:1, False:0}
RUN_EVAL=$3

####INTERNALS###################################################
################################################################
RED='\033[0;31m'
NC='\033[0m' # No Color
#scan build directory for pm-v filename format(lists all methods)
i=0
while read line
do
      binaries[ $i ]="$line"        
          (( i++ ))
done < <(find ../build/pm-v-* -printf "%f\n")


#If we are runninng only a subet of kernels,
#use the user-specified list.
if [ "$RUN_SUBSET_OF_KERNELS" -eq 1 ]; then
  echo "using subset of kernels:"
  IFS=', ' read -r -a kernellist <<< "${@:5}"
  echo $kernellist
else
  #Get all available kernels from the base binary
  KERNELS=$(../build/pm-v-gcc-O3 --listkernels d)
  IFS=', ' read -r -a kernellist <<< "$KERNELS"
fi


if [ ${#binaries[@]} -eq 0 ]; then
    printf "$RED No pm-v-* builds found in /build.$NC \n"
    printf "Use 'make' in project root to build them.\n"
  else
    # If it exists, delete timing file and create a new one.
    # All binaries write their timings into the same file
    # One column in the resulting file contains the binary name,
    # which simplifies plotting with R
    rm -f $RUN_RESULT/evaluation.txt
    touch $RUN_RESULT/evaluation.txt

    # Run evaluation of image pair within patchmatch
    if [ "$RUN_EVAL" -eq 1 ]; then
      echo "running evaluation within ./pm"
      EVALPARMS="--mask $MASK_FILE  --gt $GT_FILE"
      echo "binary_name, kernel_name, rows, cols, time, w, q, fpc, avg_mcost_time, sppp, coverage, bad05, bad1, invalid, avgErr " > $RUN_RESULT/evaluation.txt
    else #only record timing
      echo "binary_name, kernel_name, rows, cols, time, w, q, fpc, avg_mcost_time, sppp" > $RUN_RESULT/evaluation.txt
    fi

    for binversion in "${binaries[@]}"; do
      echo "Timing binary version: $binversion"

      # run each kernel that was found
      for kernel in "${kernellist[@]}"; do
        echo "Timing Kernel: $kernel"
        ADDITIONAL_PARAMS="--kernel $kernel --unit-test 0"
        BIN=../build/$binversion


        #./core/compute-single.sh $IMAGE_DIR/$IMAGE_NAME $RUN_RESULT $ADDITIONAL_PARAMS
        #Create experiment run detail file
        touch $RUN_RESULT/run-params.txt
        echo "Git commit hash:" > $RUN_RESULT/run-params.txt
        echo "$(git rev-parse HEAD)" >> $RUN_RESULT/run-params.txt 
        echo "Git branch:" > $RUN_RESULT/run-params.txt
        echo "$(git branch | grep \* | cut -d ' ' -f2)" >>  $RUN_RESULT/run-params.txt
        echo "Additional args:" >> $RUN_RESULT/run-params.txt 
        echo $ADDITIONAL_ARGS >> $RUN_RESULT/run-params.txt

        #Run Patchmatch
        COMMAND="$BIN --left $LEFT_IMG --right $RIGHT_IMG --output $RUN_RESULT $EVALPARMS $ADDITIONAL_PARAMS"
        echo $COMMAND #For debug purposes
        $COMMAND

        #Run eval
      done
    done 
fi
exit

