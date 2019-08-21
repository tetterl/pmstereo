#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

# output is to
# runs/IMAGE_NAME
# all kernels and binary version combinations write to the same directory,
# (image results will be overwritten)
# timing and eval output is appended for the entire run

#example: run all kernels on quarter resolution Teddy image, no eval
#./time-all-kernels.sh trainingQ Teddy 0 0

#example: run all only simd_v3 kernel on quarter resolution Teddy image, with eval
./time-all-kernels.sh trainingQ Teddy 1 1 simd_v2 simd_v3


