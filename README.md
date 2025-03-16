# bianchi-maass-forms
This repository is part of my PhD Thesis "Computing Bianchi-Maass Forms for Class Number 1 Non-Euclidean Rings" supervised by Solomon Friedberg at Boston College.

It computes Maass cusp forms for the Bianchi groups $`\Gamma_d = \textup{PSL}_2(\mathcal{O}_{\mathbb{Q}(\sqrt{-d})})`$ where $`d = 19, 43, 67, 163`$.
I generalized Hejhal's method to these groups, and it serves as the core method of the computation.
This repository contains C++ code for locating probable Laplace eigenvalues and providing evidence that they are genuine via a modularity check.
It also contains SageMath code for plotting the resulting forms.
I wrote these SageMath methods in order to make what I believe is the first free to use, open source, high-resolution plots of Bianchi polyhedra.
These methods build on John Cremona's `bianchi-progs` (https://github.com/JohnCremona/bianchi-progs) which contains code that implement's Swan's algorithm.

The C++ code depends on
* Eigen (https://eigen.tuxfamily.org)
* Boost (https://www.boost.org)
* archt (https://github.com/BBDJSH-LuCaNT-2023/verification-code/tree/main/phase2/archt) which relies on...
  * MPFI (https://gitlab.inria.fr/mpfi/mpfi)
  * MPFR (https://www.mpfr.org)
  * GMP (https://gmplib.org)
* Arb which is contained in FLINT (https://flintlib.org, it also depends on MPFR and GMP)
* OpenMP (probably already part of your compiler)

I used CMake as my build system. I have included my `CMakeLists.txt` as an example.

Some commentary about the directories in this repository:
* **Output** - This contains data from the Coarse, Medium, and Final searches, as well as the Final Check step. 
The subdirectory "Tested, Hand Checked, and Final" contains manually post-processed data. 
All final data included in the Thesis comes from this subdirectory.
* **Computation Logs** - This contains terminal output from when the computations were run on the Andromeda Linux Cluster at Boston College. These are extremely verbose and long. This folder has its own README.
* **Plotting** - This contains SageMath code files and a Jupyter Notebook with example code usage. There are also some animated GIFs of plots of Bianchi-Maass forms in this directory.