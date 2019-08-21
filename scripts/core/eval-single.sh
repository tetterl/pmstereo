#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.
# User-specific params(should be passed by arguments)
if [ "$#" -ne 2 ]
  then
    echo "Illegal number of arguments supplied. Usage:"
    echo "eval-single.sh IMAGE_NAME RUN_RESULT_DIR"
    exit 1
fi
# Path to MiddEval3/Teddy (or other image) directory created by download-dataset.sh
# Path to the eval binary
BIN=../../build/eval
# Where the groundtruth is stored
GT_DIR=../../dataset/MiddEval3/$1
# Path to the directory where the patch-match binary put it's output
RUN_RESULT=$2

#
PFM_FILE=$RUN_RESULT/disp1.pfm
GT_FILE=$GT_DIR/disp0GT.pfm
MASK_FILE=$GT_DIR/mask0nocc.png
EVAL_OUT=$RUN_RESULT/evaluation.txt

#directories should all exist (compute.sh created them)
$($BIN $PFM_FILE $GT_FILE $MASK_FILE > $EVAL_OUT)
echo "Wrote result to $EVAL_OUT"
