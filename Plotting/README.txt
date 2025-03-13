"fd_wireframe.sage" - Contains the function "fd_wireframe(.....)" which plots fundamental domains for Bianchi groups d = 19, 43, 67, 163 (and others may or may not work; it is not fully tested otherwise).

"KBesselPrecomputer.sage" - Contains a class which precomputes values of a K-Bessel function of order ir and then does a linear interpolation for any argument that is at least as large as the the passed in "minimum_argument". Absolute error achieved is to within something like 1e-6 at best. How many points are used per unit interval can be adjusted in the class initializer. This accuracy is definitely good enough for plotting.

"Bianchi-Maass Plots.ipynb" - Jupyter Notebook that demonstrates how to use the "fd_wireframe(....)$ function.


HOW TO RUN THESE PROGRAMS

Assuming you have SageMath and Jupyter Notebook installed:

1. Download bianchi-progs by John Cremona (https://github.com/JohnCremona/bianchi-progs)
2. Place the three files in the folder "FD" from bianchi-progs.
3. Open "Bianchi-Mass Plots.ipynb" and read the commented code which demonstrates how to use the functions.

If you want to plot Bianchi-Maass forms, you will need to access a file which contains the index and coefficient data. Methods are provided which are able to read in the data produced by my C++ program.