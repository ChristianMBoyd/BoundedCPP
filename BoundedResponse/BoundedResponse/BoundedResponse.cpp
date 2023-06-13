// BoundedResponse.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Loss.h"

int main()
{
    // surface response calculation class
    Loss g;

    // test of loss() with significant input anisotropy for comparison with Mathematica
    //std::cout << "\nThe result: ";
    //std::cout << g.loss(0.05, 0.01, 1., 10., 7., 3., 5., 1.1, 0.1, 1000, 5, 1.);

    // test of mChi0DiagNew()
    double q = 0.1, w = 1.1, delta = 0.1, Qn = 1.8, Qnp = 0.01, L = 10, cutoff = 5;
    bool evenPar = false;
    int nMax = g.nMax(cutoff, L);
    bool evenMax = g.evenQ(nMax);
    Eigen::VectorXd Qlist = g.Qlist(L, nMax, evenPar, evenMax);
    std::cout << "\n The result from mChi0Diag(): " << g.mChi0Diag(q, w, delta, Qlist, L, evenPar).diagonal();
    std::cout << "\n The result from mChi0DiagNew(): " << g.mChi0DiagNew(q, w, delta, Qlist, L, evenPar).diagonal();

    // To do:
    //  1) Double-check mChi0Diag() was implemented correctly, then move on
    //  2) Optimize mChi0OffDiag() by only calculating the non-zero entries of Pi0Qn() and Pi0Qnp()
    //      Note: no more logic checks, compare to previous cases first
    //  3) Check Eigen documentation again to optimize input/return types

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
