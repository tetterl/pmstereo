#only run on gcc-O3-time-mcost
# we need a seperate node since the level 3 cache is shared (smart cache)

bsub -n 12 -R "select[model=XeonE5_2680v3]" ./size-dependency.sh simd_v3 1 simd_v3
bsub -n 12 -R "select[model=XeonE5_2680v3]" ./size-dependency.sh simd_v17 1 simd_v17
bsub -n 12 -R "select[model=XeonE5_2680v3]" ./size-dependency.sh simd_v21 1 simd_v21

