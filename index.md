## Super-Fibonacci Spirals

[Marc Alexa](https://www.cg.tu-berlin.de/team/prof-dr-marc-alexa/)

### Abstract

Super-Fibonacci spirals are an extension of Fibonacci spirals, enabling fast generation of an arbitrary but fixed number of 3D orientations. The algorithm is simple and fast. A comprehensive evaluation comparing to other methods shows that the generated sets of orientations have low discrepancy, minimal spurious components in the power spectrum, and almost identical Voronoi volumes. This makes them useful for a variety of applications, in particular Monte Carlo sampling. 

### Documents

[CVPR paper](https://github.com/marcalexa/superfibonacci/docs/superfib_cvpr.pdf), includes supplemental material

### Code

See the repository. Dependencies: Eigen, CLI11. [Code from Julie Mitchell](https://mitchell-web.ornl.gov/SOI/index.php) is included (copyrighted  but free for commercial use).

Uses old school makefile to allow compiling sequential and parallel code with different complilers. On MacOS, clang generates faster seuqnetial code but still lacks support for C++17 standard parallelism. 

### Data

See the repository

Includes some data from Charles Karney. Complete data available on [github](https://github.com/cffk/orientation)

