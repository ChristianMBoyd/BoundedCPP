# BoundedCPP
This project contains the calculation for the surface loss function of a bounded, anisotropic electron gas ported from Mathematica to C++.

The basic linear algebra interface is handled by Eigen, which is included.  Matching Mathematica's C-compilation performance, however, requires implementing MKLsubroutines.

The MKL .dll files are excluded from this repository.  The necessary Intel oneAPI inclusions are 1) the .dll files from the oneAPI/mkl redist folder and 2) the libiomp5md.dll file from the oneAPI/compiler redist folder.



