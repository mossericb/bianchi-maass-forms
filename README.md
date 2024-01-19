# bianchi-maass-forms
Code for computing Bianchi-Maass cusp forms for class number 1 rings of integers, including the non-Euclidean rings.

The SageMath version of this computation is correct but far too slow for the sheer size of this computation if we want to compute Maass forms for fields larger than Q(i). In particular, the non-Euclidean cases all have discriminant high enough that great care must be given to code performance, not to mention it needs to be run in parallel on many cores to stand a chance.

A C++ version is written and working. Right now, it's in the stage of optimizing for memory and CPU usage via implementing a new idea for leveraging symmetry to the greatest extent possible in Hejhal's algorithm.
