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
	std::complex<double> Pi0(double q, double w, double delta, double Qn, double Qnp);
	double xiScale(double qx, double qy, double ex, double ey, double ez);
	std::complex<double> sumPi0Interval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0(double q, double w, double delta, double Qn, double L);
	int nMax(double cutoff, double L);
	bool evenQ(int val);
	Eigen::VectorXi posList(const int parity, const int nMax);
	Eigen::VectorXi intList(const int parity, const int nMax); // may be unnecessary
	Eigen::VectorXi posToIntList(const int parity, const int nMax);
	double qScale(double qx, double qy, double mx);
	double LScale(double L, double mz);
	Eigen::VectorXcd mChi0DiagList(double q, double w, double delta, Eigen::VectorXd& Qlist, double L);
	Eigen::MatrixXcd mChi0OffDiagList(double q, double w, double delta, Eigen::VectorXd& Qlist);
	bool zeroQ(double Q, double L);
	Eigen::MatrixXcd mChi0(double qx, double qy, double mx, double mz, double w, double delta, double L, double cutoff, const int parity);

	// consider putting internal functions into private
private:
	const std::complex<double> i; // sqrt[-1]
	const double pi;
};