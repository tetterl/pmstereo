#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

IMAGE_NAME=trainingQ/Teddy
RUN_RESULT=../../runs/Teddy
#To be implemented in pm: Take additional params that control
#Which subalgos are run
ADDITIONAL_PARAMS=
./core/compute-single.sh $IMAGE_NAME $RUN_RESULT $ADDITIONAL_PARAMS


