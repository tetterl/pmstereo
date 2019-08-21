#only run on gcc-O3
bsub -R "select[model=XeonE5_2680v3]" ./size-dependency.sh baseline 1 baseline 
bsub -R "select[model=XeonE5_2680v3]" ./size-dependency.sh weight_precomputation_fast_exp 1 weight_precomputation_fast_exp 
bsub -R "select[model=XeonE5_2680v3]" ./size-dependency.sh rgba 1 rgba 
bsub -R "select[model=XeonE5_2680v3]" ./size-dependency.sh simd_v7 1 simd_v7 
bsub -R "select[model=XeonE5_2680v3]" ./size-dependency.sh simd_v10 1 simd_v10
