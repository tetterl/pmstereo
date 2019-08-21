#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

IMAGE_NAME=trainingF/Teddy
RUN_RESULT=../../runs/Teddy

./core/eval-single.sh $IMAGE_NAME $RUN_RESULT
