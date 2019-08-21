#!/bin/bash
cd "${0%/*}" #change shell to script's dir. solves relative path problems.

# output is to
# runs/{trainingQ,trainingF}/IMAGE_NAME
# all kernels and binary version combinations write to the same directory,
# (image results will be overwritten)
# timing and eval output is appended for the entire run


# Use all CPUs on node to avoid other processes interfering with memory 

# Quarter resolution images, Run all kernels and all configurations. 
bsub -n 18 -W 22:00 -R "select[model=XeonE5_2697v4]" ./time-all-kernels.sh trainingQ Teddy 1 0 

