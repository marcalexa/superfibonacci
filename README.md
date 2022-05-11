## Super-Fibonacci Spirals

Fast, low discrepancy sampling of orientations.

### Code

[s3onlyfib.cc](src/s3onlyfib.cc) provides a minimal C++ implementation of Super-Fibonacci sampling without dependencies (apart from STL). It expects the number of samples as the sole command line parameter and outputs the rotations (as quaternions) to `stdout`. 

All other source code depends on [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for storing and possibly manipulating the quaternions and  [CLI11](https://github.com/CLIUtils/CLI11) for command line interpretation.

An old school makefile for compliation is provided. Unlike `cmake` this allows compiling sequential and parallel code with different complilers. On MacOS, `clang` generates faster sequential code but still lacks support for C++17 standard parallelism. Makefile requires modification to find compilers, dependencies, as usual. 

All tools write and read orientations as quaternions, in ASCII of double precision floating point numbers. 

Slightly modified code based on  [Julie Mitchell's implementation](https://mitchell-web.ornl.gov/SOI/index.php) of the SOI and Hopf fibration methods is included (copyrighted  but free for commercial use).

### Data

Data directory contains various pre-generated rotation samples and shell scripts to generate statistics provided in the paper. 

Includes some data from Charles Karney. Complete data available on [github](https://github.com/cffk/orientation)

