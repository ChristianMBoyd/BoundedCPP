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
    //std::cout << "\nThe result of posToIntList: ";
    //for (int item : g.posToIntList(1, 13)) // nice C++ 11 feature
    //{
    //    std::cout << item << ", ";
    //}

    //std::cout << "\nThe result of intList: ";
    //for (int item : g.intList(1, 13)) // nice C++ 11 feature
    //{
    //    std::cout << item << ", ";
    //}

    //std::cout << "\nThe result of posList(posToIntList) ";
    //for (int item : g.posToIntList(1, 13)) // nice C++ 11 feature
    //{
    //    Eigen::VectorXi posList = g.posList(1, 13);
    //    std::cout << posList[item] << ", ";
    //}

    // tests of Loss::mChi0OffDiagList()
   /* Eigen::VectorXd Qlist = (3.14 / 5) * g.posList(0, g.nMax(5, 5)).cast<double>();
    Eigen::MatrixXcd offDiagList = g.mChi0OffDiagList(0.1, 1.1, 0.1, Qlist);
    std::cout << "\nThe initialized part of offDiagList: " << offDiagList(2,1);
    std::cout << "\n The uninitialized part of offDiagList: " << offDiagList(1,2);*/

    // test of Loss::mChi0 -- OLD VERSION!  Currently new inputs, check defs before using!
    /*Eigen::MatrixXcd mChi0 = g.mChi0(0.1, 0.0, 1., 1., 1.1, 0.1, 4, 1000, 0); // DEPRECATED INPUT
    std::cout << "\nThe result of mChi0 at specified position: " << mChi0(199,200);*/

    // test of (-1) powers operating as expected
   /* std::cout << "\nThe result of (-1)^0 is: " << pow(-1, 0);
    std::cout << "\nThe result of (-1)^1 is: " << pow(-1, 1);*/

    // test of mChi0New
    //const int parity = 1;
    //const double L = 4.0;
    //const int nMax = g.nMax(5, L); // cutoff = 5
    //Eigen::VectorXd Qlist = (3.141592653589793238463 / L) * g.posList(parity, nMax).cast<double>();
    //std::cout << "\nThe result under investigation:\n";
    //Eigen::MatrixXcd mChi0 = g.mChi0New(0.1, 1.1, 0.1, Qlist, L, parity);
    //Eigen::MatrixXcd mChi0OffDiag = g.mChi0OffDiag(0.1, 1.1, 0.1, Qlist);
    //Eigen::MatrixXcd mChi0Diag = g.mChi0Diag(0.1, 1.1, 0.1, Qlist, L, parity);
    //Eigen::MatrixXcd test(2, 2);
    //test << 1, 2, 3, 4;

    // test of loss(), using mChi0New() -- NEXT: ImChi --- find ways to reduce Eigen:: initializations!
    std::cout << "The result: ";
    std::cout << g.loss(0.1, 0.0, 1., 1., 1., 1., 1., 1.1, 0.1, 50, 5, 1.);

    return 0;
}
