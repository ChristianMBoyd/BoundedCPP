// BoundedResponse.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include "Loss.h"

int main()
{
    // surface response calculation class
    Loss g;

    // test of loss() with significant input anisotropy for comparison with Mathematica
    std::cout << "\nThe result: ";
    std::cout << g.loss(0.001, 0.001, 1., 100., 1., 3., 5., 0.1, 0.1, 500, 5, 1.);

    // parametric tests of different mChi0 implementations --- seem to agree to double precision
    /*double q = 0.001, w = 0.1, delta = 0.1, L = 500, cutoff = 5;
    const bool evenPar = true;
    const int nMax = g.nMax(cutoff, L);
    const bool evenMax = g.evenQ(nMax);
    Eigen::VectorXd Qlist = g.Qlist(L, nMax, evenPar, evenMax);
    Eigen::MatrixXcd oldM = g.mChi0DiagOld(q, w, delta, Qlist,L,evenPar);
    Eigen::MatrixXcd newM = g.mChi0Diag(q, w, delta, Qlist, L, evenPar);
    Eigen::MatrixXcd diff = (oldM - newM).diagonal();
    std::cout << "The sum of differences is:\n";
    std::cout << diff.sum();*/

    // To do:
    //  1) Slight difference with Mathematica at small q and large mass/dielectric anisotropy
    //  2) Check Eigen documentation again to optimize input/return types
    //  3) Go through and correct *OR* try to understand better "auto" usage
    //      i.e., likely remove auto calls on Eigen objects

    //  Next:
    //  1) Before CUDA, do Mathematica tests on ParallelTable vs. Table
    //      i.e., is it better to parallelize the spectrum or the calculation
    //  2) Make a second project for CUDA implementation for backend

    // closing preamble, allows output code to hang for comparison/inspection in release
    char input;
    std::cout << "\nEnter any input to close.\n";
    std::cin >> input;

    return 0;
}
