# add list of all kernels to run behind simd_v3
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_baseline 1 baseline
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_split_up_cost 1 split_up_cost 
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_weight_precomputation_fast_exp 1 weight_precomputation_fast_exp 
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_rgba 1 rgba 
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_simd_v3 1 simd_v3 
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_simd_v7 1 simd_v7 
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_simd_v10 1 simd_v10
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-kernel.sh mcost_simd_v15 1 simd_v15
