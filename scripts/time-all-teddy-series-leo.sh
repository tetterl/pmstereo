# add list of all kernels to run behind simd_v3
bsub -R "select[model=XeonE5_2680v3]" ./time-teddy-series.sh split_up_cost 1 split_up_cost 
bsub -R "select[model=XeonE5_2697v4]" ./time-teddy-series.sh weight_precomputation_fast_exp 1 weight_precomputation_fast_exp 
bsub -R "select[model=XeonE5_2697v4]" ./time-teddy-series.sh rgba 1 rgba 
bsub -R "select[model=XeonE5_2697v4]" ./time-teddy-series.sh simd_v3 1 simd_v3 
bsub -R "select[model=XeonE5_2697v4]" ./time-teddy-series.sh simd_v7 1 simd_v7 
bsub -R "select[model=XeonE5_2697v4]" ./time-teddy-series.sh simd_v10 1 simd_v10
