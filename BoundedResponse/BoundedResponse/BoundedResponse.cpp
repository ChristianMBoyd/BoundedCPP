// BoundedResponse.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Loss.h"

int main()
{
    // surface response calculation class
    Loss g;

    // test of loss() with significant input anisotropy for comparison with Mathematica
    std::cout << "\nThe result: ";
    std::cout << g.loss(0.05, 0.01, 1., 10., 7., 3., 5., 1.1, 0.1, 1000, 5, 1.);

    // Note: manual omp calls seem to have had little effect on the release build, but debug and debugging the release are quicker
    //      i.e., likely that the compiler has already optimized the for loops

    // To do:
    //  1) Enumerate only the non-zero entries BEFORE calling Pi0
    //  2) Then, pull |Qn|,|Qnp|<1 logic out of Pi0 and call specific parts
    //      i.e., one matrix is built from the Qn argument, the other from Qnp, then they're combined
    //  3) Make a second project for CUDA implementation for backend

    // closing preamble, allows output code to hang for comparison/inspection in release
    char input;
    std::cout << "\nEnter any input to close.\n";
    std::cin >> input;

    return 0;
}
