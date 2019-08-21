#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.
#We require 3 parameters (if running all kernels)
#else more
if [ "$#" -lt 1 ]
  then
    echo "Illegal number of arguments supplied. Usage:"
    echo "./time-teddy-kernel.sh JOB_NAME RUN_SUBSET kernels..."
    echo "output will be put in dir /runs/TeddySeries/JOB_NAME"
    echo "RUN_SUBSET = 1 -> run subset of kernels specified"
    echo "RUN_SUBSET = 0 -> run all kernels"
    echo "example:"
    echo "./time-teddy-series.sh 0 "
    exit 1
fi

RUN_RESULT=../runs/TeddySeries/$1

# Create output directory if it does not exist
mkdir -p $RUN_RESULT

#Eval paths and variables
EVAL_BIN=../build/eval

#Directory with the im0, im1 ground truth images
SCALESERIES_DIR=../dataset/TeddyScales

#If this is set to 0, ALL kernels are run!
RUN_SUBSET_OF_KERNELS=$2 #{True:1, False:0}

####INTERNALS###################################################
################################################################
RED='\033[0;31m'
NC='\033[0m' # No Color
#scan build directory for pm-v filename format(lists all methods)
i=0
while read line
do
  echo $line
  binaries[ $i ]="$line"        
  (( i++ ))
done < <(find ../build/pm-v-* -printf "%f\n")


#If we are runninng only a subet of kernels,
#use the user-specified list.
if [ "$RUN_SUBSET_OF_KERNELS" -eq 1 ]; then
  echo "using subset of kernels:"
  IFS=', ' read -r -a kernellist <<< "${@:3}"
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
    echo "binary_name, kernel_name, rows, cols, time, w, q, fpc, sppp" > $RUN_RESULT/evaluation.txt

    for binversion in "${binaries[@]}"; do
      echo "Timing binary version: $binversion"

      # run each kernel that was found
      for kernel in "${kernellist[@]}"; do
        echo "Timing Kernel: $kernel"
        ADDITIONAL_PARAMS="--kernel $kernel --unit-test 0"
        BIN=../build/$binversion
        for SCALE in 100
        do 
          LEFT_IMG=$SCALESERIES_DIR/im0-$SCALE.png
          RIGHT_IMG=$SCALESERIES_DIR/im1-$SCALE.png
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
          $COMMAND | tee $RUN_RESULT/"$(basename $BIN)"_std_out.txt #actually run command

        #Run eval
        done 
      done
    done 
fi
exit

