# bianchi-maass-forms
This program is part of my PhD thesis work.

It computes Bianchi-Maass cusp forms for class number 1 rings of integers, including the non-Euclidean rings (the fields Q(sqrt(-d)) for d = 1, 2, 3, 7, 11, 19, 43, 67, 163).

Please excuse the cobwebs; I'm in the stage of converging on the final version of the program and dumping large commented blocks, unused methods, etc. 

This contains methods for locating probable eigenvalues, running heuristic checks on the results, evaluating, and plotting the resulting Bianchi-Maass forms. Functionality for verifying various conjectures such as Ramanujan-Petersson and Sato-Tate is coming soon.

This program relies on the libraries Eigen, OpenMP, and OpenCV (for plot generation, can be safely removed otherwise). It also relies on archt from https://github.com/BBDJSH-LuCaNT-2023/verification-code.git (and its dependencies mpfi, mpfr, gmp). While archt computes rigorous bounding intervals for K-Bessel values, my program does not use interval arithmetic, but the archt methods are sufficiently fast and come with an accuracy guarantee which is still useful for operating in double precision.

The SageMath version of this computation is correct and functional but far too slow to compute Maass forms for fields larger than the first few. In particular, the non-Euclidean cases all have discriminant high enough that a highly-parallel C++ version is required.