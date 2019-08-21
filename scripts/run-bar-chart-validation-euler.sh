# Run this on gcc-O3 binary only!
bsub -R "select[model=XeonE5_2680v3]" ./time-all-kernels-final.sh trainingQ Teddy 1 0 

