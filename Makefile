PHONY: all

all: allgcc allintel

clean:
	rm build/pm-v-*

#ALLINTEL:

allgcc: GCC-O3 GCC-03-ftvec GCC-03-funsafe GCC-03-ffast-math GCC-O3-time-mcost
allintel:  ICC-Ofast ICC-O3 

GCC-O3: pm.cpp
	g++ -std=c++11 -O3 -march=native -mavx2 -mfma -o build/pm-v-gcc-O3 pm.cpp -lpng
GCC-O2: pm.cpp
	g++ -std=c++11 -O2 -march=native -mavx2 -mfma -o build/pm-v-gcc-O2 pm.cpp -lpng
GCC-O0: pm.cpp
	g++ -std=c++11 -O0 -march=native -mavx2 -mfma -o build/pm-v-gcc-O0 pm.cpp -lpng
# for timing the mcost function, adds a little overhead
GCC-O3-time-mcost: pm.cpp
	g++ -DTIME_MCOST -std=c++11 -O3 -march=native -mavx2 -mfma -o build/pm-v-gcc-O3-time-mcost pm.cpp -lpng


GCC-03-ftvec: pm.cpp
	g++ -std=c++11 -O3 -march=native -mavx2 -mfma -ftree-vectorize -o build/pm-v-gcc-O3-ftvec pm.cpp -lpng
GCC-03-funsafe: pm.cpp
	g++ -std=c++11 -O3 -march=native -mavx2 -mfma -funsafe-math-optimizations -o build/pm-v-gcc-O3-funsafe pm.cpp -lpng
GCC-03-ffast-math: pm.cpp
	g++ -std=c++11 -O3 -march=native -mavx2 -mfma -ffast-math -o build/pm-v-gcc-O3-ffast-math pm.cpp -lpng

#ICC
#fma is enabled by the -xavx2  flag automatically
ICC-Ofast: pm.cpp
	icc -std=c++11 -Ofast -axAVX,CORE-AVX2 -o build/pm-v-iccOfast pm.cpp -lpng
ICC-O3: pm.cpp
	icc -std=c++11 -O3 -axAVX,CORE-AVX2 -o build/pm-v-icc-O3 pm.cpp -lpng
ICC-O2: pm.cpp
	icc -std=c++11 -O2 -axAVX,CORE-AVX2 -o build/pm-v-icc-O2 pm.cpp -lpng
ICC-O: pm.cpp
	icc -std=c++11 -O  -axAVX,CORE-AVX2 -o build/pm-v-icc-O pm.cpp -lpng



