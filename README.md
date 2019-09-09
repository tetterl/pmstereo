# Patch match stereo implementation

Implements the Patch Match stereo algorithm by Michael Bleyer, Christoph Rhemann and Carsten Rother presented in [PatchMatch Stereo - Stereo Matching with Slanted Support Windows](https://microsoft.com/en-us/research/wp-content/uploads/2011/01/PatchMatchStereo_BMVC2011_6MB.pdf). Implemented and optimized at ETH 2019 during the [How to Write Fast Numerical Code Lecture](https://acl.inf.ethz.ch/teaching/fastcode/2019/) by Nik Bamert, Jordan Burklund and Thomas Etterlin.

## Abstract

High-performance implementation of PatchMatch Stereo matching algorithm on the Intel Haswell architecture. Demonstrating an 8.9x speedup over the C++ baseline implementation by using methods such as arithmetic transformations, precomputations, access pattern optimization, SIMD usage, etc.
Roofline analysis and an instruction latency analysis shows that the implementation is within a factor of 2 of the optimal performance.

## Report

See https://github.com/tetterl/pmstereo/blob/master/report.pdf

## Attribution 
Note that the very first baseline implementation was based on the following repository: https://github.com/ivanbergonzani/patch-match-stereo.

For the baseline implementation, the entire code from Ivan Bergonzani's was rewritten to remove the OpenCV library dependency.
Our baseline implementation uses only standard C++ and STL features and the PNG library to load and store images.
Please note that some variable names from the original repository may remain but that the code has been extensively reworked.

Further note that for evaluation purposes, the Middlebury evaluation toolkit has been partially integrated into the codebase. The original
source code can be obtained from http://vision.middlebury.edu/stereo/code/. 


# Setup
To run pmstereo, some stereo ground truth images are required. Please acquire them by running the following commands,
which download and extract the dataset into the subdirectory ```/dataset```.

```
/scripts/core/download-dataset-fullres.sh 
/scripts/core/download-dataset-half.sh 
/scripts/core/download-dataset.sh 
```

## Build instructions (default binary, cmake build system)
The following builds a default binary without a specific compiler or flags in mind.

```
cd pmstereo
mkdir build
cd build
cmake ..
make
build/pmstereo ...flags... 
```
For a list of flags, please see subsection flags below.

## Build instructions (gcc,icc, flags)
The following builds binaries with various compiler flags for evaluation by the scripts within the subdirectory /script.
The targets are all build within the /build subdirectory and all binary names are prefixed with the string ```pm-v-````.
The intention behind this naming scheme is that the scripts in /scripts can find all binary versions by themselves without
passing an excplicit binary list.

```
cd pmstereo
make

```

# Usage 

## Flags 
All flags are prefixed with --, i.e. 

```
./pmstereo --left image1.png --right image2.png 
```

| Flag        | Values          | Description                                                               |
|-------------|-----------------|---------------------------------------------------------------------------|
| left        | file  path      | path to left input image                                                  |
| right       | file path       | path to right input image                                                 |
| output      | directory path  | path to where the output (disparity, logs) of the run should be written.  |
| gt          | file path       | ground truth image (only if evaluation included in run)                   |
| mask        | file path       | mask ground truth image (only if evaluation included in run)              |
| kernel      | string          | name of the kernel to be run                                              |
| listkernels | string          | lists available kernels (left, right path still have to be supplied)      |
| maxdisp     | integer         | the maximum disparity to use for the matching                             |
| unit-test   | integer (0,1)   | whether or not to run the unit test of the cost function on the kernels   |
| sppp        | integer         | single precision floating point peak performance of current machine       |


