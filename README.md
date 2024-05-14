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
| Signal 1 | 0.113338s  | 0.00527319s | 0.956901s		  | 0.0094s	 | 0.982385s	   | 0.00204517s  | 9.64566s              | 0.111637s  |
| Signal 2 | 0.104097s | 0.00335322s | 0.96229s		  | 0.0112685s  | 0.975212s		  | 0.014643s  | 9.62468s      | 0.14903s |
| Signal 3 | 0.103579s | 0.00145639s | 0.961407s		  | 0.0107727s  | 0.968442s	  | 0.0114715s | 9.64459s      | 0.164847s |

The test.cpp file can be used to compare results with the PyWavelets library. The results are similar, but not exact. They should be more accurate as this algorithm
uses the integral to calculate the exact value at a point, rather than an integral estimation method used in PyWavelets.

# Usage

This library currently supports 1D vectors, although I do have plans to create a 2D version. Downloading the `.cpp` and `.h` files should be enough to run in your environment.
Refer to `test.cpp` for example usage. 

# Dependencies

I say dependencies, although the code used is basic enough to the point that it may work on older versions. The dependencies listed are those that 
I've used in creating this library.

- CMake (>=3.21)
- C++ 17 Compiler
- OpenMP (>= 5)

# Benchmark Reproduction

The benchmark file can be compiled using the CMakeLists.txt file. To reproduce benchmark results, the CWT.exe file can be run with argument values for the number of signal samples and the number of frequencies 
(i.e: `$ CWT.exe 10000 300`). 