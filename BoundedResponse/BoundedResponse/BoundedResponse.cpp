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
    std::cout << g.loss(0.05, 0.01, 1., 10., 7., 3., 5., 1.1, 0.1, 500, 5, 1.);

    // To do:
    //  implement CUDA

    // closing preamble, allows output code to hang for comparison/inspection in release
    char input;
    std::cout << "\nEnter any input to close.\n";
    std::cin >> input;

    return 0;
}
