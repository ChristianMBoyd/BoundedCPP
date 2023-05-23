// BoundedResponse.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Loss.h"

int main()
{
    Loss g;
    double x, y;
    std::cout << "Hello World!";

    // test of new branch cut in Loss::Sqrt()
    std::cout << "\nPlease enter the real part: ";
    std::cin >> x;
    std::cout << "\nPlease enter the imaginary part: ";
    std::cin >> y;
    std::complex<double> z(x, y);
    std::cout << "\nThe value entered was: " << z;
    std::cout << "\nThe standard complex square root of the value entered is: " << std::sqrt(z);
    std::cout << "\nThe adjusted complex square root of the value entered is: " << g.posRoot(z);

    // test of Loss::Pi0Int() accuracy
    std::cout << "\nThe result of Pi0Int: " << g.Pi0Int(0.4, 0.1, 0.1, 0.3, 0.5);

    // test of Loss::Xi() accuracy
    std::cout << "\nThe result of Xi: " << g.Xi(0.1, 0.2, 1, 3, 5);

    // test of Loss::sumPi0() accuracy -- needs to be fixed
    std::cout << "\nThe result of sumPi0: " << g.sumPi0(0.1, 1.1, 0.1, 0.4, 100);

    return 0;
}

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
