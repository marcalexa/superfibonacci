## Super-Fibonacci Spirals

Fast, low discrepancy sampling of orientations.

### Code

Dependencies: [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [CLI11](https://github.com/CLIUtils/CLI11). Slightly modified code based on  [Julie Mitchell's implementation](https://mitchell-web.ornl.gov/SOI/index.php) of the SOI and Hopf fibration methods is included (copyrighted  but free for commercial use).

Uses old school makefile to allow compiling sequential and parallel code with different complilers. On MacOS, clang generates faster sequential code but still lacks support for C++17 standard parallelism. Makefile requires modification to find compilers, dependencies, as usual. 

All tools write and read orientations as quaternions, in ASCII of double precision floating point numbers. 


### Data

Includes some data from Charles Karney. Complete data available on [github](https://github.com/cffk/orientation)

Data directory also contains shell scripts to generate data provided in the paper. 
