# bianchi-maass-forms
This program is part of my PhD thesis work.

It computes Bianchi-Maass cusp forms for class number 1 rings of integers, including the non-Euclidean rings (the fields Q(sqrt(-d)) for d = 1, 2, 3, 7, 11, 19, 43, 67, 163).

This contains methods for locating probable eigenvalues, running heuristic checks on the results, evaluating, and plotting the resulting Bianchi-Maass forms. Functionality for verifying various conjectures such as Ramanujan-Petersson and Sato-Tate is coming soon.

Searching to find Bianchi-Maass cusp forms is a large computation that needs to be run across many cores to complete in a reasonable amount of time. 

This program relies on many libraries, including Arb (and its dependencies mpfi, mpfr, Flint, gmp), Eigen, OpenMP, and OpenCV (for plot generation, can be safely removed otherwise).

The SageMath version of this computation is correct and functional but far too slow to compute Maass forms for fields larger than the first few. In particular, the non-Euclidean cases all have discriminant high enough that a highly-parallel C++ version is required.
