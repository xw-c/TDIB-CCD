# TDIB-CCD

Code to run the [Wang21 benchmark](https://archive.nyu.edu/handle/2451/61518). 

## Included Libraries
This repo contains the following libraries:
- [CCD-Query-IO](https://github.com/Continuous-Collision-Detection/CCD-Query-IO): Handles input/output operations for the benchmark.
- [rational-cpp](https://github.com/zfergus/rational-cpp): Wrapper of GMP's Rational Type.
- [GMP](https://gmplib.org/): Reading rational numbers from strings.
  - A precompiled dynamic library of GMP 6.3.0 is provided for Windows users. The source code is in [here](https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz).
  - For other platforms (Linux, macOS), you need to install GMP manually. Refer to the GMP official documentation for installation instructions.

## Usage
### Compilation
To compile the code, simply run:
```
xmake b benchmark
```

### Running the Benchmark
Run the benchmark on wang21 dataset with:
```
xmake r benchmark -p <path-to-benchmark>
```
Replace `<path-to-benchmark>` with the path to your Wang21 benchmark dataset.

## Notes
- The precompiled GMP library in this repository is only for Windows. If you're using a different platform, please install GMP manually.
- `xmake` manages all other dependencies, so no additional installations are required for the included libraries.
