// BoundedResponse.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Loss.h"

int main()
{
    Loss g;
    /*double x, y;
    std::cout << "Hello World!";*/

    // test of new branch cut in Loss::Sqrt()
   /* std::cout << "\nPlease enter the real part: ";
    std::cin >> x;
    std::cout << "\nPlease enter the imaginary part: ";
    std::cin >> y;
    std::complex<double> z(x, y);
    std::cout << "\nThe value entered was: " << z;
    std::cout << "\nThe standard complex square root of the value entered is: " << std::sqrt(z);
    std::cout << "\nThe adjusted complex square root of the value entered is: " << g.posRoot(z);*/

    // test of Loss::Pi0Int() accuracy
   /* std::cout << "\nThe result of Pi0Int: " << g.Pi0Int(0.4, 0.1, 0.1, 0.3, 0.5);*/

    // test of Loss::Xi() accuracy
   /* std::cout << "\nThe result of Xi: " << g.Xi(0.1, 0.2, 1, 3, 5);*/

    // test of Loss::sumPi0() accuracy -- needs to be fixed
    //std::cout << "\nThe result of sumPi0: " << g.sumPi0(0.1, 1.1, 0.1, 0.4, 100);
    
    // test of Loss::posTointList()  --- built with posList()
    std::cout << "\nThe result of posToIntList: ";
    for (int item : g.posToIntList(1, 13)) // nice C++ 11 feature
    {
        std::cout << item << ", ";
    }

    std::cout << "\nThe result of intList: ";
    for (int item : g.intList(1, 13)) // nice C++ 11 feature
    {
        std::cout << item << ", ";
    }

    std::cout << "\nThe result of posList(posToIntList) ";
    for (int item : g.posToIntList(1, 13)) // nice C++ 11 feature
    {
        Eigen::VectorXi posList = g.posList(1, 13);
        std::cout << posList[item] << ", ";
    }

    return 0;
}
