# BemTool
BemTool is a C++ header library that only relies on the STL and Boost (www.boost.org/).
Its main purpose is the computation of assembly routines for boundary elements. At low level,
it contains a set of elementary classes for low dimensional linear algebra (so as to deal with
geometric computations) and simplicial meshes in 1D,2D and 3D. It contains basic wrappers
for Gmsh (cf http://gmsh.info/) as regards mesh generation, and eigen (http://eigen.tuxfamily.org/)
for the treatment of large matrices.

Currently it provides routines for the assembly boundary element matrices for Laplace,
Helmholtz, Yukawa (purely dissipative Hekmholtz) and Maxwell matrices in 2D and 3D.
The available discretisation spaces are P0, continuous P1, continuous P2 and lowest order
Raviart-Thomas elements.

One of the motivations for the development of this lilbrary was numerical investigation
of the performance of the multi-trace formalism. As a consequence the classes for meshes
include functionalities specially oriented toward the management of multi-domain geometries,
in particular for cases where subdomains are adjacent so as to give rise to junction points.

#### The library content
Besides a makefile and .cpp files at the top level of the library, BemTool
is composed of three directories. The main directory is bemtool/ which contains 
all the (header only) source code. The other two directories (doc/ and mesh/) are
not strictly necessary for the library to run. The directory mesh/ contains a series
of toy meshes generated with Gmsh. The directory doc/ contains .tex and .pdf files
that describes the diagonalisation of classical boundary integral operators by means
of separation of variables on a unit circle (in 2D) or sphere (in 3D).

#### Calling the library
In BemTool, a special attention has been paid to the ease of use. Linking with the library
only requires a \#include instruction toward the file tools.hpp located in the directory
bemtool of the library.

#### Exemple files
There is currently no documentation available, but test2D.cpp and test3D.cpp are exemple files
that show how to call the library. These should be compilable through commands make test2D
and make 3D.