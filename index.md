## Super-Fibonacci Spirals

[Marc Alexa](https://www.cg.tu-berlin.de/team/prof-dr-marc-alexa/)

### Abstract

Super-Fibonacci spirals are an extension of Fibonacci spirals, enabling fast generation of an arbitrary but fixed number of 3D orientations. The algorithm is simple and fast. A comprehensive evaluation comparing to other methods shows that the generated sets of orientations have low discrepancy, minimal spurious components in the power spectrum, and almost identical Voronoi volumes. This makes them useful for a variety of applications, in particular Monte Carlo sampling. 

### Documents

[CVPR paper](https://github.com/marcalexa/superfibonacci/files/8650030/superfib_cvpr.pdf), includes supplemental material

### Erratum

Super Fibonacci Spirals, as described in the paper, have no simple refinement property: The paper claims that a set of kn samples generated with Super-Fibonacci sampling contains the set generated for n samples and, more generally, sets of m and n samples would have every k-th sample in common, where k is the gcd of m and n. This is not true. 

### Code

The main point of Super-Fibonacci sampling for generating rotations is that it is as simple as uniform random sampling but on par with more sophisticated methods in terms of the properties of the distribution. In pseudo-code, the procedure for generating *n* rotations represented as quaternions *q<sub>i</sub>* is:

> For *i* in (0,...,*n*-1)
> 1. *s* = *i* + 0.5
> 2. r = sqrt(*s*/*n*)
> 3. R = sqrt(1.0-*s*/*n*)
> 4. &alpha; = (2&pi; * *s*)/&phi;
> 5. &beta; = (2&pi; * *s*)/&psi;
> 6. q<sub>i</sub> = (r sin α, r cos α, R sin β, R cos β)

Here, &phi; and &psi; are two magic constants, for which I recommend the values &phi; = &radic;2 and &psi; = 1.533751168755204288118041...

This pseudocode should be trivial to convert into actual code. Examples in python

```
import numpy as np

phi = np.sqrt(2.0)
psi = 1.533751168755204288118041

Q = np.empty(shape=(n,4), dtype=float)

for i in range(n):
    s = i+0.5
    r = np.sqrt(s/n)
    R = np.sqrt(1.0-s/n)
    alpha = 2.0 * np.pi * s / phi
    beta = 2.0 * np.pi * s / psi
    Q[i,0] = r*np.sin(alpha)
    Q[i,1] = r*np.cos(alpha)
    Q[i,2] = R*np.sin(beta)
    Q[i,3] = R*np.cos(beta)
```

and C/C++ (including some minor optimizations for speed)

```
std::vector<std::array<double,4> > Q(n);

double dn = 1.0 / (double)n;
double mc0 = 1.0 / sqrt(2.0);
double mc1 = 1.0 / 1.533751168755204288118041;

for (int i = 0; i < n; i++)
{
    s = (double)i+0.5;
    ab = 2.0 * M_PI * s;
    alpha = ab * mc0;
    beta = ab * mc1;
    s *= dn;
    r = sqrt(s);
    R = sqrt(1.0-s);
    Q[i][0] = r*sin(alpha);
    Q[i][1] = r*cos(alpha);
    Q[i][2] = R*sin(beta);
    Q[i][3] = R*cos(beta);
}
```

The repository contains more code, allowing to experiment with magic constants, genrate other distributions, and provides tools for analysis of the sampling. 


### Data

See the repository for precomputed sequences of various sampling strategies and shell scripts to generate more data or compare the sample sets. 

