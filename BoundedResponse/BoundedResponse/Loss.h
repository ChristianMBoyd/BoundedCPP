#pragma once
#include <cmath>
#include <complex>
#include <vector> // need variable-size arrays for below
#include "Eigen/Dense" // need linear algebra tools

class Loss
{
public: 
	Loss();
	std::complex<double> posRoot(std::complex<double>);
	std::complex<double> Pi0Int(double q, double w, double delta, double Qn, double Qnp);
	double Xi(double qx, double qy, double ex, double ey, double ez);
	std::complex<double> sumPi0Interval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0(double q, double w, double delta, double Qn, double L);
	int nMax(double cutoff, double L);
	bool evenQ(int val);
	Eigen::VectorXi posList(const int parity, const int nMax);
	Eigen::VectorXi intList(const int parity, const int nMax); // may be unnecessary
	Eigen::VectorXi posToIntList(const int parity, const int nMax); 
	Eigen::MatrixXcd mChi0(double qx, double qy, double mx, double mz, double w, double delta, double L, double cutoff, const int parity);

private:
	const std::complex<double> i; // sqrt[-1]
	const double pi;
};