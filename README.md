# Ricker Continuous Wavelet Transform (RCWT)

Create this C++ library as there seems to be little to no CWT libraries
for the Ricker Wavelet. Note that this library is poor relative to many other available libraries.
The benchmark times can be found below. The current version is somewhat basic as I've translated code
from the PYWT library. I've made a few optimisations within this code, namely shortening the convolution to
only what's required, removing the usage of the `j` array used in PYWT. Using multi-threading with OpenMP, I've been able to reduce
calculation times significantly. Going forward, I think a move is required towards a Fast Fourier Transform, but again, 
this is not something I'm familiar with so I'll need to figure that out first.

# Benchmark

Tests have been carried out on an ASUS Zenbook with 8 GB 2400MHz SODIMM, 2.1 GHz AMD Ryzen 5 3500U with Radeon Vega Mobile Gfx.

## Results

Below are the benchmark results for each signal. This can be carried out independently using the benchmark.cpp file.

| Signal   | 10k - 300  | SD 		       | 10k - 3000 	 | SD 		        | 100k - 300 	 | SD 		       | 100k - 3000 	 | SD 		     |
|----------|------------|-------------|--------------|--------------|--------------|-------------|---------------|-----------|
| Signal 1 | 0.267584s  | 0.00403663s | 2.54619s		  | 0.0101705s	 | 2.57582s	   | 0.0175132s  | 27.2342s      | 0.201316s  |
| Signal 2 | 0.272648s | 0.00832729s | 2.54661s		  | 0.00406283s  | 2.58572s		  | 0.0157627s  | 27.3185s      | 0.18166s |
| Signal 3 | 0.267189s | 0.00316097s | 2.54882s		  | 0.0119903s  | 2.58616s		  | 0.0134953s | 27.2408s      | 0.14322s |

The test.cpp file can be used to compare results with the PyWavelets library. The results are similar, but not exact. If you require more precision,
consider changing the sample size to a larger value. The default is 1000 samples of the Ricker Wavelet across the bounds [-8, 8]. 

# Usage

This library currently supports 1D vectors, although I do have plans to create a 2D version.

# Dependencies

I say dependencies, although the code used is basic enough to the point that it may work on older versions. The dependencies listed are those that 
I've used in creating this library.

- CMake (>=3.21)
- C++ 17 Compiler
- OpenMP (>= 5)

# Benchmark Reproduction

The benchmark file can be compiled using the CMakeLists.txt file. To reproduce benchmark results, the CWT.exe file can be run with argument values for the number of signal samples and the number of frequencies 
(i.e: `$ CWT.exe 10000 300`). 